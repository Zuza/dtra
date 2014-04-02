#include <cctype>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <mutex>
#include <queue>
#include <utility>
#include <vector>
#include "core/mapping.h"
#include "core/klcs.h"
#include "core/lis.h"
#include "core/util.h"
#include "core/Coverage.h"
#include "ssw/ssw_cpp.h"

#include "core/BoundedStringDistance.h"
using namespace std;

// ovdje saram s intovima i size_t-ovima, iako je
// 32 bitni int dovoljan svuda

DEFINE_string(long_read_algorithm, "coverage", "Algorithm for long reads (naive|lis|coverage|ssw)");
DEFINE_bool(multiple_hits, true, "Allow multiple placements on a single gene.");
DEFINE_double(windowed_alignment_size, 2, "Windowed alignment size.");

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
                         const int seedLen,
                         const shared_ptr<Read>& read,
                         const int rc) {
    int noN = 0; // sluzi da bih mogao preskociti seedove s N-ovima, Y-onima i sl.
    
    typedef unsigned long long ullint;
    thread_local vector<pair<ullint, ullint> > positions;
    assert(seedLen < 256-1);
    thread_local char read_sub[256];

    for (int i = 0; i < read->size(); ++i) {
      if (!isBase(read->get(i, rc))) { ++noN; }
      if (i >= seedLen && !isBase(read->get(i-seedLen, rc))) { --noN; }
    
      if (noN == 0 && i+1 >= seedLen) {
        for (int j = 0; j < seedLen; ++j)
          read_sub[j] = read->get(i-seedLen+1+j, rc);

        idx->get_substring_pos(positions, read_sub, seedLen);
        for (auto pos_pair : positions) {
          pair<unsigned int, unsigned int> x(pos_pair.first, pos_pair.second);
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
                         const int seedLen,
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
               positions[start].first + FLAGS_windowed_alignment_size * read->size()
               > positions[end].first) {
          ++end;
        }

        assert(start <= end);

        // TODO(fpavetic): After klcs stabilize, remove
        //                 lis and coverage.
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
          vector<Interval> intervals;
          for (int i = start; i < end; ++i) {
            int a = positions[i].first;
            int b = positions[i].first + seedLen - 1;
            int c = positions[i].second;
            intervals.push_back(Interval(a,b,c));
          }

          if (!intervals.empty()) {
            int score = 0;
            vector<Interval> coverage;
            cover(&coverage, &score, intervals);
	  
            assert(!coverage.empty());
            int begin = coverage[0].left - coverage[0].value;
            int end = coverage.back().right+read->size()-coverage.back().value;
	  
            read->updateMapping(score, begin, end, rc, geneIdx,
                                genes[geneIdx]->description(), "");	
          }
        } else if (FLAGS_long_read_algorithm == "klcs") {
          vector<pair<int, int> > match_pairs(positions.begin()+start,
                                              positions.begin()+end);
          if (!match_pairs.empty()) {
            int gene_offset = match_pairs[0].first;
            int mx1 = 0, mx2 = 0;
            for (size_t i = 0; i < match_pairs.size(); ++i) {
              match_pairs[i].first -= gene_offset;
              assert(match_pairs[i].first >= 0);
              mx1 = max(mx1, match_pairs[i].first);
              mx2 = max(mx2, match_pairs[i].second);
            }

            int klcs_score = 0;
            vector<pair<int, int> > klcs_recon;
            klcs(match_pairs, seedLen, &klcs_score, &klcs_recon);

            int beg = gene_offset + klcs_recon[0].first - klcs_recon[0].second;
            // TODO: za sada se ne koristi, ali treba testirati izraz za end
            int end = gene_offset + klcs_recon.back().first - 
              read->size() - klcs_recon.back().second;

            read->updateMapping(klcs_score, beg, end, rc, geneIdx,
                                genes[geneIdx]->description(), "");
          }
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
      populatePositions(positionsByGene, genes, idx, idx->getSeedLen(), read, rc);
    
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
        windowedAlignment(read, positionsByGene, idx->getSeedLen(), genes, rc);
      }
    }
  }
  
}

template<int MAXD>
void calcEditDistance(vector<shared_ptr<Gene> >& genes,
                      shared_ptr<Read> read) {
  // TODO: trenutno se ovaj objekt kreira za svaki read, idealno bi bilo
  //       kad bi si ga dretva radilica napravila jednom i samo racunala
  //       readove koji je dopadnu
  BoundedStringDistance<false, MAXD> bsd(1);

  assert(!genes.empty());

  for (auto it = read->topMappings().begin(); it != read->topMappings().end(); ++it) {
    int geneIdx = it->geneIdx;
    string description = it->geneDescriptor;
    int editDistance = -1;
    
    if (geneIdx < genes.size() &&
        genes[geneIdx]->description() == description) { // nuzna provjera jer
                                                        // geneIdx je indeks na
                                                        // razini bloka, a ne
                                                        // cijele baze
      // kandidat je u trenutno procesuiranom bloku!

      const char* text = genes[geneIdx]->data() + it->geneBegin;
      int patternLen = read->size();

      int textWindowSize = it->geneEnd-it->geneBegin+1;
      double score = it->score; // score je mjera identicnih poklapanja
      
      // gornja granica na gresku je kad i u textu i u patternu ostavimo
      // identicna poklapanja i obrisemo sve ostalo
      int limit = (int)(fabs(textWindowSize-score)+fabs(patternLen-score));
      limit = min(limit, MAXD-1);
      limit = max(limit, 1);
      
      char* pattern = new char[patternLen+1];
      read->toCharArray(pattern, it->isRC);

      editDistance = bsd.compute(text, pattern, patternLen, limit);

      delete[] pattern;
      it->editDistance = editDistance;
    }
  }
}

// * ako je fillEditDistance false, pokrece se 'pametni' algoritam procjena
// regija u kojima se nalazi read
// * ako je fillEditDistance true, podrazumijeva se da je read vec popunjen
// procjenama regija i treba samo izracunati i dopuniti informaciju o 
// edit-distanceu
// * nekad je dovoljno samo procijeniti regiju, dok je pun rezim rada zamisljen
// kao poziv ove funkcije s fillEditDistance = false nad svim blokovima fasta
// baze na ulazu programa. Zatim se ponovno nad svim blokovima baze poziva ova
// funkcija s fillEditDistance = true
void performMapping(vector<shared_ptr<Gene> >& genes,
                    shared_ptr<Index> idx, 
                    shared_ptr<Read> read,
                    bool fillEditDistance) {
  if (!fillEditDistance) {
    performMappingLong(genes, idx, read);
  } else {
    // najgori slucaj podrazumijevam 10% greske pa tako podesavam MAXD
    if (read->size() <= 100) {
      calcEditDistance<10>(genes, read);
    } else if (read->size() <= 500) {
      calcEditDistance<50>(genes, read);
    } else if (read->size() <= 1000) {
      calcEditDistance<100>(genes, read);
    } else if (read->size() <= 2000) {
      calcEditDistance<200>(genes, read);
    } else {
      calcEditDistance<300>(genes, read);
    }
  }
}
