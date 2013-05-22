#include <cctype>
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

// readovi duljine manje od ove obradjuju se algoritmom za krace readovce,
// ostali algoritmom za dulje
const int kShortLongBorder = 0;

DEFINE_string(long_read_algorithm, "coverage", "Algorithm for long reads (lis|coverage)");

namespace {

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

    // a) prvo procijenim koje su pozicije najizglednije
    // da na njih smjestam read

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
    //    preko odredjenog postotka hitova
    vector<pair<int, pair<int, int> > > revBegEst;
    for (auto iter : beginEstimate) {
      revBegEst.push_back(make_pair(iter.second, iter.first));
    }
    sort(revBegEst.rbegin(), revBegEst.rend());
    
    const double usableFraction = 0.6;
    int currCount = 0;
    for (size_t i = 0; i < revBegEst.size(); ++i) {
      currCount += revBegEst[i].first;
      if (1.0 * currCount / totalCount > usableFraction) {
    	while (revBegEst.size() > i+1) { revBegEst.pop_back(); }
    	break;
      }
    }

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
        if (lastStart + 1.7 * read->size() > start) {
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
}

}

void performMapping(vector<shared_ptr<Gene> >& genes,
		    shared_ptr<Index> idx, 
		    shared_ptr<Read> read) {
  assert(FLAGS_long_read_algorithm == "lis" || 
	 FLAGS_long_read_algorithm == "coverage");

  if (read->size() < kShortLongBorder) {

  } else {
    performMappingLong(genes, idx, read);
  }
}
