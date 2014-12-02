#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>
#include <memory>

#include "core/gene.h"
#include "core/util.h"
#include "core/klcs.h"
#include "core/lis.h"


using namespace std;

static void printUsageAndExit() {
  printf("glis \n");
  exit(1);
}

static vector<shared_ptr<Gene> > readGenesFromFile(FILE* inputFile) {
  assert(inputFile);

  vector<shared_ptr<Gene> > genes;
  while (true) {
    shared_ptr<Gene> g(new Gene());
    
    if (!readGene(g.get(), inputFile)) {
      break;
    } else {
      genes.push_back(g);
    }
  }

  // for (size_t i = 0; i < genes.size(); ++i) {
  //   printGene(genes[i].get(), stdout);
  // }

  return genes;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    // fasta1 fasta2
    fprintf(stderr, "Wrong format: ./glis <fasta1> <fasta2>\n");
    exit(1);
  }

  vector<shared_ptr<Gene> > A = readGenesFromFile(fopen(argv[1], "r"));
  vector<shared_ptr<Gene> > B = readGenesFromFile(fopen(argv[2], "r"));

  printf("Input reading finished... %lu %lu\n", A[0]->dataSize(), B[0]->dataSize());

  A[0]->subNonAcgtWithRandom();
  B[0]->subNonAcgtWithRandom();

  printf("Starting forward strand...\n");
  {
    vector<pair<int, int> > matches;
    int klcs_len; klcs(A[0]->data(), A[0]->dataSize(), B[0]->data(), B[0]->dataSize(), 31, &klcs_len, &matches);
    printf("matches.size() = %lu\n", matches.size());

    vector<int> result;
    calcLongestIncreasingSubsequence(&result, matches);

    printf("result.size() = %lu\n", result.size());
    FILE* out = fopen("fwd.out", "w");
    for (int idx : result) {
      fprintf(out, "%d %d\n", matches[idx].first, matches[idx].second);
    }
    fclose(out);
  }

  printf("Forward strand finished. Starting reverse...\n");
  {
    A[0]->reverseAndComplement();

    vector<pair<int, int> > matches;
    int klcs_len; klcs(A[0]->data(), A[0]->dataSize(), B[0]->data(), B[0]->dataSize(), 31, &klcs_len, &matches);
    printf("matches.size() = %lu\n", matches.size());

    vector<int> result;
    calcLongestIncreasingSubsequence(&result, matches);

    printf("result.size() = %lu\n", result.size());
    FILE* out = fopen("rev.out", "w");
    for (int idx : result) {
      fprintf(out, "%d %d\n", (int)A[0]->dataSize() - matches[idx].first, matches[idx].second);
    }
    fclose(out);
  }

  printf("Finished\n");
  return 0;
}
