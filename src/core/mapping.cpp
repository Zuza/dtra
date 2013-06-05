#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <mutex>
#include <queue>
#include <utility>
#include <vector>
#include "core/mapping.h"
#include "core/lis.h"
#include "core/util.h"
#include "core/Coverage.h"
using namespace std;

// ovdje saram s intovima i size_t-ovima, iako je
// 32 bitni int dovoljan svuda

DEFINE_string(long_read_algorithm, "lis", "Algorithm for long reads (naive|lis|coverage)");
DEFINE_bool(multiple_hits, false, "Allow multiple placements on a single gene.");

const int kBeginEstimateGroupDist = 15;

namespace {

void printUsageAndExit() {
 printf("For single hits: --long_read_algorithm=lis|naive\n");
 printf("For multiple hits: --multiple_hits --long_read_algorithm=lis|coverage\n");
 exit(1);
}

int estimateBeginFromLis(const vector<pair<int, int> >& positions,
			 const vector<int>& lis) {
  map<int, int> beginEstimate;
  for (size_t i = 0; i < lis.size(); ++i) {
    int r = lis[i];
    int p = positions[r].first;
    int kmer = positions[r].second;

    ++beginEstimate[max(0, p-kmer)];
  }
  
  int maxBegin = 0;
  int begin = -1;
  for (auto it : beginEstimate) {
    if (it.second > maxBegin) {
      maxBegin = it.second;
      begin = it.first;
    }
  }
  return begin;
}

void singleLis(shared_ptr<Read> read, 
	       const map<int, shared_ptr<vector<pair<int, int> > > >& positionsByGene,
	       const vector<shared_ptr<Gene> >& genes,
	       const int rc) {
  for (auto candidateGenes : positionsByGene) {
    int geneIdx = candidateGenes.first;
    
    vector<pair<int, int> >& positions = *candidateGenes.second;
    sort(positions.begin(), positions.end());
    
    vector<int> lisResult;
    calcLongestIncreasingSubsequence(&lisResult, positions);
    
    // gdje procjenjujemo da je pocetna pozicija
    // mapiranja reada na gen?
    int begin = estimateBeginFromLis(positions, lisResult);
      
    // segment na genu gdje procjenjujemo mapiranje
#ifdef DEBUG
    string geneSegment = cstrToString(genes[geneIdx]->data() + begin,
				      read->size());
#else
    string geneSegment = "";
#endif
    
    read->updateMapping(lisResult.size(), begin, rc, geneIdx,
			genes[geneIdx]->description(), geneSegment);
  }
}

void createBeginEstimates(
  map<int, shared_ptr<map<int, int> > >& beginEstimateByGene,
  const map<int, shared_ptr<vector<pair<int, int> > > >& positionsByGene
) {
  for (auto candidateGenes : positionsByGene) {
    int geneIdx = candidateGenes.first;

    vector<pair<int, int> >& positions = *candidateGenes.second;
    shared_ptr<map<int, int> > beginEstimate(new map<int, int>());

    for (size_t i = 0; i < positions.size(); ++i) {
      int position = positions[i].first;
      int kmer = positions[i].second;
      
      if (position-kmer >= 0) {
	++((*beginEstimate)[position-kmer]);
      }
    }

    beginEstimateByGene[geneIdx] = beginEstimate;
  }
}

void groupNeighbouringEstimates(
    map<int, shared_ptr<map<int, int> > >& beginEstimateByGeneGrouped,
    const map<int, shared_ptr<map<int, int> > >& beginEstimateByGene) {
  for (auto be : beginEstimateByGene) {
    int geneIdx = be.first;
    shared_ptr<map<int, int> >& beginEstimate = be.second;
    shared_ptr<map<int, int> > beginEstimateGrouped(new map<int, int>());

    for (map<int, int>::iterator it = beginEstimate->begin(); 
	 it != beginEstimate->end(); ++it) {
      int total = 0;
      
      for (map<int, int>::iterator prev = it; 
	   prev != beginEstimate->begin(); ) {
	--prev;
	if (it->first-prev->first <= kBeginEstimateGroupDist) {
	  total += it->second;
	} else {
	  break;
	}
      }
     
      map<int, int>::iterator next = it; ++next;
      for ( ; next != beginEstimate->end(); ++next) {
	if (next->first-it->first <= kBeginEstimateGroupDist) {
	  total += it->second;
	} else {
	  break;
	}
      }

      (*beginEstimateGrouped)[it->first] = total;
    }

    beginEstimateByGeneGrouped[geneIdx] = beginEstimateGrouped;
  }
}

void singleNaive(shared_ptr<Read> read, 
		 const map<int, shared_ptr<vector<pair<int, int> > > >& positionsByGene,
		 const vector<shared_ptr<Gene> >& genes,
		 const int rc) {
  map<int, shared_ptr<map<int, int> > > beginEstimateByGene;
  createBeginEstimates(beginEstimateByGene, positionsByGene);
  
  map<int, shared_ptr<map<int, int> > > beginEstimateByGeneGrouped;
  // groupNeighbouringEstimates(beginEstimateByGeneGrouped,
  // 			       beginEstimateByGene);
    
  beginEstimateByGeneGrouped = beginEstimateByGene;
    
  // TODO: trenutno ovo vrijedi, nece nakon updatea grupiranja
  assert(beginEstimateByGene.size() == beginEstimateByGeneGrouped.size());

  if (beginEstimateByGene.empty()) {
    assert(positionsByGene.empty());
    return;
  }

  // naive -> procjena pozicije naivnim brojanjem
  for (auto iter : beginEstimateByGeneGrouped) {
    int geneIdx = iter.first;
    shared_ptr<map<int, int> > beginEstimateGrouped = iter.second;
    
    int bestHits = 0;
    int bestPos = -1;
    for (map<int, int>::iterator it = beginEstimateGrouped->begin();
	 it != beginEstimateGrouped->end(); ++it) {
      if (it->second > bestHits) {
	bestPos = it->first;
	bestHits = it->second;
      }
    }
    
    assert(bestHits > 0);
    read->updateMapping(bestHits, bestPos, rc, geneIdx,
			genes[geneIdx]->description(), "");
  }
}

void populatePositions(map<int, shared_ptr<vector<pair<int, int> > > >& positionsByGene,
		       const vector<shared_ptr<Gene> >& genes,
		       const shared_ptr<Index>& idx,
		       const shared_ptr<Read>& read,
		       const int rc) {
  unsigned long long hsh = 0;
  int seedLen = idx->getSeedLen();
  unsigned long long andMask = (1LL<<(2*seedLen))-1;
  int noN = 0; // sluzi da bih mogao preskociti seedove s N-ovima, Y-onima i sl.
  
  for (int i = 0; i < read->size(); ++i) {
    if (!isBase(read->get(i, rc))) { ++noN; }
    if (i >= seedLen && !isBase(read->get(i-seedLen, rc))) { --noN; }
    
    hsh = hsh*4+baseToInt(read->get(i, rc));
    hsh &= andMask;
    
    if (noN == 0 && i+1 >= seedLen) {
      Index::iterator it = idx->getPositions(hsh, seedLen);
      for ( ; !it.done(); it.advance()) {
	pair<unsigned int, unsigned int> x = it.get();
	int geneId = x.first;
	int position = x.second;
	//printf("geneId=%d, position=%d\n", geneId, position);
	
	shared_ptr<vector<pair<int, int> > >& posVec =
	  positionsByGene[geneId];
	
	if (!posVec) {
	  posVec = shared_ptr<vector<pair<int, int> > > (
	  	      new vector<pair<int, int> >());
	}

	assert(positionsByGene.count(geneId));
	posVec->push_back(make_pair(position, i+1-seedLen));
      } 
    }
  }
}

// TODO: ova se funkcija moze srediti tako da koristi 
// helper funkcije napisane u ovom fileu,
// trenutno je ovo samo kopi pejst nekog starog commita
// PROBLEM: cini se da ne radi trenutno bas
void windowedAlignment(shared_ptr<Read> read,
		       map<int, shared_ptr<vector<pair<int, int> > > >& positionsByGene,
		       const shared_ptr<Index> idx,
		       const vector<shared_ptr<Gene> >& genes,
		       const int rc) {
  // a) prvo procijenim koje su pozicije najizglednije
  // da na njih smjestam read
  int seedLen = idx->getSeedLen();

  map<pair<int, int>, int> beginEstimate;
  int totalCount = 0;

  for (auto candidateGenes : positionsByGene) {
    int geneIdx = candidateGenes.first;
    
    vector<pair<int, int> >& positions = *candidateGenes.second;
    for (size_t i = 0; i < positions.size(); ++i) {
      int position = positions[i].first;
      int kmer = positions[i].second;
      const int blockSize = 20;
      int block = (position-kmer)/blockSize;
      ++beginEstimate[make_pair(geneIdx, max(0,block*blockSize))];
      ++totalCount;
    }
  }
  
  // b) odredim prozore za provjeru, a to su oni koji imaju zajedno
  // preko odredjenog postotka hitova
  vector<pair<int, pair<int, int> > > revBegEst;
  for (auto iter : beginEstimate) {
    revBegEst.push_back(make_pair(iter.second, iter.first));
  }
  sort(revBegEst.rbegin(), revBegEst.rend());
  
  map<int, shared_ptr<vector<int> > > candidatePositionsByGene;
  for (size_t i = 0; i < revBegEst.size(); ++i) {
    int geneIdx = revBegEst[i].second.first;
    int pos = revBegEst[i].second.second;
    shared_ptr<vector<int> >& vecPtr =
      candidatePositionsByGene[geneIdx];
    if (!vecPtr) {
      vecPtr = shared_ptr<vector<int> > (new vector<int>());
    }
    vecPtr->push_back(pos);
  }
  
  for (auto candidatePositions : candidatePositionsByGene) {
    int geneIdx = candidatePositions.first;
    vector<int>& starts = *candidatePositions.second;
    sort(starts.begin(), starts.end());
    
    vector<pair<int, int> >& positions = *positionsByGene[geneIdx];
    sort(positions.begin(), positions.end());
    
    int windowSize = 2*read->size();
    int b = 0, e = 0;
    int lastStart = -1000000000;

    for (auto start : starts) {
      if (lastStart + read->size() > start) {
	continue;
      }
      lastStart=start;
      
      while (b < positions.size() &&
	     positions[b].first < start) { ++b; }
      while (e < positions.size() &&
	     positions[e].first < start+windowSize) { ++e; }
      assert(b <= e);

      int score = 0;
      if (FLAGS_long_read_algorithm == "lis") {
        vector<pair<int, int> > pripremaZaLis;
	for (int i = b; i < e; ++i) {
	  pripremaZaLis.push_back(positions[i]);
	}

	vector<int> lisResult;
        calcLongestIncreasingSubsequence(&lisResult, pripremaZaLis);
	score = lisResult.size();
      } else if (FLAGS_long_read_algorithm == "coverage") {
	vector<Interval> intervals;
	for (int i = b; i < e; ++i) {
	  int a = positions[i].first;
	  int b = positions[i].first + seedLen - 1;
	  int c = positions[i].second;
	  intervals.push_back(Interval(a,b,c));
	}
	cover(NULL, &score, intervals);
      }

#ifdef DEBUG
      string geneSegment = cstrToString(genes[geneIdx]->data() + begin,
					read->size());
#else
      string geneSegment = "";
#endif

      read->updateMapping(score, start, rc, geneIdx,
			  genes[geneIdx]->description(), geneSegment);

    }
  }
}

void performMappingLong(vector<shared_ptr<Gene> >& genes,
			shared_ptr<Index> idx, shared_ptr<Read> read) {
  for (int rc = 0; rc < 2; ++rc) {
    map<int, shared_ptr<vector<pair<int, int> > > > positionsByGene;
    populatePositions(positionsByGene, genes, idx, read, rc);
    
    if (!FLAGS_multiple_hits) {
      if (FLAGS_long_read_algorithm == "lis") {
	singleLis(read, positionsByGene, genes, rc);
      } else if (FLAGS_long_read_algorithm == "naive") {
	singleNaive(read, positionsByGene, genes, rc);
      } else {
	printUsageAndExit();
      }
    } else { // --multiple_hits
      // TODO: ovo je trenutno potrgano, ako stignem popravljam
      // windowedAlignment(read, positionsByGene, idx, genes, rc);
    }
  }
}
  
}

void performMapping(vector<shared_ptr<Gene> >& genes,
		    shared_ptr<Index> idx, 
		    shared_ptr<Read> read) {
  performMappingLong(genes, idx, read);
}
