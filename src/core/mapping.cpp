#include <cctype>
#include <algorithm>
#include <map>
#include <mutex>
#include <utility>
#include <vector>
#include "core/mapping.h"
#include "core/lis.h"
#include "core/util.h"

using namespace std;

// ovdje saram s intovima i size_t-ovima, iako je
// 32 bitni int dovoljan svuda

// readovi duljine manje od ove obradjuju se algoritmom za krace readovce,
// ostali algoritmom za dulje
const int kShortLongBorder = 0;

namespace {

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

}

void performMapping(vector<shared_ptr<Gene> >& genes,
		    shared_ptr<Index> idx, shared_ptr<Read> read) {
  if (read->size() < kShortLongBorder) {

  } else {
    performMappingLong(genes, idx, read);
  }
}
