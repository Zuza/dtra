#include <cctype>
#include <algorithm>
#include <map>
#include <mutex>
#include <utility>
#include <vector>
#include "core/mapping.h"
#include "core/lis.h"
#include "core/util.h"

#include <gflags/gflags.h>
#include "ssw/ssw_cpp.h"
#include <iostream>

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
        vector<pair<unsigned int, unsigned int> > positions;
	
	// TODO: mozda u Index mogu staviti metodu koja direktno puni
	// positionsByGene kako bih izbjegao kopiranje podataka
	
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

DEFINE_int32(ssw_match, 1, "ssw match");
DEFINE_int32(ssw_mismatch, 2, "ssw mismatch");
DEFINE_int32(ssw_gap_open, 5, "ssw gap open");
DEFINE_int32(ssw_gap_extend, 2, "ssw gap extend");

void performSswMapping(vector<shared_ptr<Gene> >& genes,
			shared_ptr<Index> idx, shared_ptr<Read> read) {

  using namespace StripedSmithWaterman;

  shared_ptr<Read> read_rc(new Read());
  *read_rc = *read;
  read_rc->complement();

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

    aligner.Align(read->data().c_str(),
                  genes[geneIdx]->data(), genes[geneIdx]->dataSize(),
                  filter, &alignment);

    read->updateMapping(alignment.sw_score,
                        alignment.ref_begin, /*rc*/ 0, geneIdx,
                        genes[geneIdx]->description(), "");

    aligner.Align(read_rc->data().c_str(),
                  genes[geneIdx]->data(), genes[geneIdx]->dataSize(),
                  filter, &alignment);

    read->updateMapping(alignment.sw_score,
                        alignment.ref_begin, /*rc*/ 1, geneIdx,
                        genes[geneIdx]->description(), "");
  }

}

}

void performMapping(vector<shared_ptr<Gene> >& genes,
		    shared_ptr<Index> idx, shared_ptr<Read> read) {

  performSswMapping(genes, idx, read);

  // if (read->size() < kShortLongBorder) {

  // } else {
  //   performMappingLong(genes, idx, read);
  // }
}
