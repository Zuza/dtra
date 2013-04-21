#include <cctype>
#include <algorithm>
#include <map>
#include <mutex>
#include <utility>
#include <vector>
#include "core/mapping.h"
#include "core/lis.h"
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
      if (!isBase(toupper(read->get(i, rc)))) { ++noN; }
      if (i >= seedLen && !isBase(toupper(read->get(i-seedLen, rc)))) { --noN; }

      hsh = hsh*4+baseToInt(toupper(read->get(i, rc)));
      hsh &= andMask;
      
      if (noN == 0 && i+1 >= seedLen) {
        vector<pair<unsigned int, unsigned int> > positions;
	
        // vector kojeg vrati ova metoda je sortiran poretkom
        // u kojem se pairovi standardno sortiraju
        idx->getPositions(&positions, hsh);
	
        for (auto x : positions) {
          int geneId = x.first;
          int position = x.second;
          
          // TODO: positions vector je sortiran kao sto se pairovi inace
          // sortiraju pa mozda mogu izbjeci trazenje po mapi svaki put
          if (!positionsByGene.count(geneId)) {
            positionsByGene[geneId] = 
              shared_ptr<vector<pair<int, int> > > (new vector<pair<int, int> >);
          }
          positionsByGene[geneId]->push_back(make_pair(position, i+1-seedLen));
        }
      }
    }

    // ovo se ne bi smjelo dogoditi,
    // barem na sintetski generiranim podacima
    // assert(!positionsByGene.empty());
    
    for (auto candidateGenes : positionsByGene) {
      int geneId = candidateGenes.first;
      vector<pair<int, int> >& positions = *candidateGenes.second;

      // treba li sort
      sort(positions.begin(), positions.end());

      vector<int> lisResult;
      calcLongestIncreasingSubsequence(&lisResult, positions);
      
      int begin = calcBegin(positions, lisResult);
      // printf("geneId=%d, begin=%d, score=%d\n", 
      // 	     geneId, begin, (int)lisResult.size());
      string geneSegment = cstrToString(genes[geneId]->data() + begin, read->size());
      string geneName = string(genes[geneId]->name() + 1);
      read->updateMapping(lisResult.size(), geneId, begin, rc, 
                          geneName, geneSegment);
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
