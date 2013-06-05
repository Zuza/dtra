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
DEFINE_bool(single_hit, false, "Assume that read maps to only one segment of a gene");

const int kBeginEstimateGroupDist = 15;

namespace {

void printUsageAndExit() {
 printf("For single hits: --single_hit --long_read_algorithm=lis|naive\n");
 printf("For multiple hits: --long_read_algorithm=lis|naive|coverage\n");
 exit(1);
}

int calcBegin(const vector<pair<int, int> >& positions,
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
  
void performMappingOld(vector<shared_ptr<Gene> >& genes,
		       shared_ptr<Index> idx, shared_ptr<Read> read) {
  unsigned long long hsh = 0;
  int seedLen = idx->getSeedLen();
  unsigned long long andMask = (1LL<<(2*seedLen))-1;
  
  for (int rc = 0; rc < 2; ++rc) {
    map<int, shared_ptr<vector<pair<int, int> > > > positionsByGene;
    int noN = 0; // sluzi da bih mogao preskociti seedove s N-ovima,
    // Y-onima i sl.
    
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
    printf("%d\n", (int)positionsByGene.size());
    
    for (auto candidateGenes : positionsByGene) {
      int geneIdx = candidateGenes.first;
      
      vector<pair<int, int> >& positions = *candidateGenes.second;
      sort(positions.begin(), positions.end());
      
      vector<int> lisResult;
      calcLongestIncreasingSubsequence(&lisResult, positions);
      
      // gdje procjenjujemo da je pocetna pozicija
      // mapiranja reada na gen?
      int begin = calcBegin(positions, lisResult);
      
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


void performMappingLong(vector<shared_ptr<Gene> >& genes,
			shared_ptr<Index> idx, shared_ptr<Read> read) {
  unsigned long long hsh = 0;
  int seedLen = idx->getSeedLen();
  unsigned long long andMask = (1LL<<(2*seedLen))-1;
 
  for (int rc = 0; rc < 2; ++rc) {
    map<int, shared_ptr<vector<pair<int, int> > > > positionsByGene;
    int noN = 0; // sluzi da bih mogao preskociti seedove s N-ovima,
                 // Y-onima i sl.
        
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

    printf("%d\n", (int)positionsByGene.size());

    map<int, shared_ptr<map<int, int> > > beginEstimateByGene;
    createBeginEstimates(beginEstimateByGene,
			 positionsByGene);

    map<int, shared_ptr<map<int, int> > > beginEstimateByGeneGrouped;
    groupNeighbouringEstimates(beginEstimateByGeneGrouped,
			       beginEstimateByGene);

    // printf("%d %d %d\n", 
    // 	   (int)positionsByGene.size(),
    // 	   (int)beginEstimateByGene.size(),
    // 	   (int)beginEstimateByGeneGrouped.size());
    if (beginEstimateByGene.empty()) {
      assert(positionsByGene.empty());
      return;
    }

    // naive -> procjena pozicije naivnim brojanjem
    if (FLAGS_long_read_algorithm == "naive") {
      for (auto iter : beginEstimateByGeneGrouped) {
	int geneIdx = iter.first;
	shared_ptr<map<int, int> > beginEstimateGrouped = iter.second;
	
	if (FLAGS_single_hit) {
	  int bestHits = 0;
	  int bestPos = -1;
	  for (map<int, int>::iterator it = beginEstimateGrouped->begin();
	       it != beginEstimateGrouped->end(); ++it) {
	    if (it->second > bestHits) {
	      bestPos = it->first;
	      bestHits = it->second;
	    }
	  }
	  read->updateMapping(bestHits, bestPos, rc, geneIdx,
			      genes[geneIdx]->description(), "");
	} else {
	  printf("not yet supported\n");
	  exit(1);
	}
      }
    } else {

    }
  }
}

}

void performMapping(vector<shared_ptr<Gene> >& genes,
		    shared_ptr<Index> idx, 
		    shared_ptr<Read> read) {
  if (FLAGS_single_hit) {
    if (FLAGS_long_read_algorithm != "lis" &&
	FLAGS_long_read_algorithm != "naive") {
      printUsageAndExit();
    }
  } else { // multiple hits
    if (FLAGS_long_read_algorithm != "lis" &&
	FLAGS_long_read_algorithm != "coverage" &&
	FLAGS_long_read_algorithm != "naive") {
      printUsageAndExit();
    }
  }

  //performMappingOld(genes, idx, read);
  performMappingLong(genes, idx, read);
}
