#include <cctype>
#include <utility>
#include <vector>
#include "core/mapping.h"
#include "core/lis.h"
using namespace std;

// readovi duljine manje od ove obradjuju se algoritmom za krace readovce,
// ostali algoritmom za dulje
const int kShortLongBorder = 0;

namespace {

void performMappingLong(MappingResult* mappingResult,
			const Database& db, shared_ptr<Read> read) {
  unsigned long long hsh = 0;
  unsigned long long andMask = (1LL<<(2*db.getSeedLen()))-1;

  int seedLen = db.getSeedLen();

  vector<vector<size_t> > positions(db.getNoGenes());
  vector<vector<size_t> > lisResults(db.getNoGenes());

  for (int rc = 0; rc < 1; ++rc) {
    int noN = 0;
    for (int i = 0; i < read->size(); ++i) {
      if (toupper((*read)[i]) == 'N') {
	++noN;
      }
      if (i >= seedLen && toupper((*read)[i-seedLen]) == 'N') {
	--noN;
      }

      hsh = hsh*4+baseToInt((*read)[i]);
      hsh &= andMask;
      
      if (noN == 0 && i+1 >= db.getSeedLen()) {
	vector<pair<int, int> > positions;
	
	// TODO: kakve su mi garancije na poredak rezultata, mozda cu morati
	// sortirati?
	// TODO: ubaciti ovo; db.getPositions(&positions, hsh);
	for (auto x : positions) {
	  int geneId = x->first;
	  int position = x->second;
	  positions[geneId].push_back(position);
	}
      }
    }

    for (int i = 0; i < db.getNoGenes(); ++i) {
      calcLongestIncreasingSubsequence(&lisResults[i], positions[i]);
    }
  }
}

}

void performMapping(MappingResult* mappingResult,
		    const Database& db, shared_ptr<Read> read) {
  if (read->size() < kShortLongBorder) {

  } else {
    performMappingLong(mappingResult, db, read);
  }
}
