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
#include "ssw/ssw_cpp.h"
using namespace std;

// ovdje saram s intovima i size_t-ovima, iako je
// 32 bitni int dovoljan svuda

DEFINE_string(long_read_algorithm, "lis", "Algorithm for long reads (naive|lis|coverage|ssw)");
DEFINE_bool(multiple_hits, false, "Allow multiple placements on a single gene.");

const int kBeginEstimateGroupDist = 15;

namespace {

void printUsageAndExit() {
 printf("For single hits: --long_read_algorithm=lis|naive|ssw\n");
 printf("For multiple hits: --multiple_hits --long_read_algorithm=lis|coverage\n");
 exit(1);
}

// TODO: malo popraviti, trenutno je implementirana dosta jednostavna
// verzija, moguce je napraviti da bolje handle-a indele
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

// funkcija izmjenjuje vektore sadrzane u positionsByGene,
// oni ce postati sortirani
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
    read->updateMapping(lisResult.size(), begin, begin+read->size(),
			rc, geneIdx,
			genes[geneIdx]->description(), "");
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

void singleNaive(shared_ptr<Read> read, 
		 const map<int, shared_ptr<vector<pair<int, int> > > >& positionsByGene,
		 const vector<shared_ptr<Gene> >& genes,
		 const int rc) {
  map<int, shared_ptr<map<int, int> > > beginEstimateByGene;
  createBeginEstimates(beginEstimateByGene, positionsByGene);
  
  if (beginEstimateByGene.empty()) {
    assert(positionsByGene.empty());
    return;
  }

  // naive -> procjena pozicije naivnim brojanjem
  for (auto iter : beginEstimateByGene) {
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
    read->updateMapping(bestHits, bestPos, bestPos+read->size(),
			rc, geneIdx,
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

void windowedAlignment(shared_ptr<Read> read,
		       map<int, shared_ptr<vector<pair<int, int> > > >& positionsByGene,
		       const shared_ptr<Index> idx,
		       const vector<shared_ptr<Gene> >& genes,
		       const int rc) {
  for (auto candidateGenes : positionsByGene) {
    int geneIdx = candidateGenes.first;
    
    vector<pair<int, int> >& positions = *candidateGenes.second;
    sort(positions.begin(), positions.end());

    int start = 0, end = 0, lastStart = -1000000000;
    
    while (true) {
      while (start < positions.size() && 
	     lastStart + read->size() > positions[start].first) {
	++start;
      }

      if (start < positions.size()) {
	lastStart = positions[start].first;
      } else {
	break;
      }

      while (end < positions.size() && 
	     positions[start].first + 2 * read->size() > positions[end].first) {
	++end;
      }

      assert(start <= end);

      if (FLAGS_long_read_algorithm == "lis") {
	vector<pair<int, int> > lisPrepare;
	for (int i = start; i < end; ++i) {
	  lisPrepare.push_back(positions[i]);
	}

	vector<int> lisResult;
	calcLongestIncreasingSubsequence(&lisResult, lisPrepare);

	// gdje procjenjujemo da je pocetna pozicija
	// mapiranja reada na gen?
	int begin = estimateBeginFromLis(lisPrepare, lisResult);
	read->updateMapping(lisResult.size(), begin, begin+read->size(),
			    rc, geneIdx, genes[geneIdx]->description(), "");
      } else if (FLAGS_long_read_algorithm == "coverage") {
	int seedLen = idx->getSeedLen();
	vector<Interval> intervals;
	for (int i = start; i < end; ++i) {
	  int a = positions[i].first;
	  int b = positions[i].first + seedLen - 1;
	  int c = positions[i].second;
	  intervals.push_back(Interval(a,b,c));
	}

	int score = 0;
	vector<Interval> coverage;
	cover(&coverage, &score, intervals);

	assert(!coverage.empty());
	int begin = coverage[0].left - coverage[0].value;
	int end = coverage.back().right+read->size()-coverage.back().value;

	read->updateMapping(score, begin, end, rc, geneIdx,
			    genes[geneIdx]->description(), "");	
      }
    }
  }
}

DEFINE_int32(ssw_match, 1, "ssw match");
DEFINE_int32(ssw_mismatch, 2, "ssw mismatch");
DEFINE_int32(ssw_gap_open, 5, "ssw gap open");
DEFINE_int32(ssw_gap_extend, 2, "ssw gap extend");

void singleSsw(shared_ptr<Read> read, 
	       const map<int, shared_ptr<vector<pair<int, int> > > >& positionsByGene,
	       const vector<shared_ptr<Gene> >& genes,
	       const int rc) {

  shared_ptr<Read> read_orig(new Read());
  *read_orig = *read;
  read_orig->removeAllLower();

  for (size_t geneIdx = 0; geneIdx < genes.size(); ++geneIdx) {
    // ----------- SSW part -------------------
    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner(FLAGS_ssw_match,
                                          FLAGS_ssw_mismatch, 
                                          FLAGS_ssw_gap_open,
                                          FLAGS_ssw_gap_extend);
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    // Aligns the query to the ref

    aligner.Align(read_orig->data().c_str(),
                  genes[geneIdx]->data(), genes[geneIdx]->dataSize(),
                  filter, &alignment);

    read->updateMapping(alignment.sw_score,
                        alignment.ref_begin, alignment.ref_end,
			rc, geneIdx,
                        genes[geneIdx]->description(), "");
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
      } else if (FLAGS_long_read_algorithm == "ssw") {        
        singleSsw(read, positionsByGene, genes, rc);
      } else {
        printUsageAndExit();
      }
    } else { // --multiple_hits
      windowedAlignment(read, positionsByGene, idx, genes, rc);
    }
  }
}
  
}

void performMapping(vector<shared_ptr<Gene> >& genes,
		    shared_ptr<Index> idx, 
		    shared_ptr<Read> read) {
  performMappingLong(genes, idx, read);
}
