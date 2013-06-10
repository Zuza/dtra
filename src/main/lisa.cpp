#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <vector>
#include <memory>

#include "core/gene.h"
#include "core/util.h"
using namespace std;

void printUsageAndExit() {
  printf("lisa dbshuffle inputdb.fasta outputdb.fasta\n");
  exit(1);
}

void shuffleDatabase(const char* inputDb,
		     const char* outputDb) {
  FILE* inputFile = fopen(inputDb, "rt");
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

  fclose(inputFile);

  srand(1603);
  random_shuffle(genes.begin(), genes.end());

  FILE* outputFile = fopen(outputDb, "wt");
  assert(outputFile);

  for (size_t i = 0; i < genes.size(); ++i) {
    printGene(genes[i].get(), outputFile);
  }
  
  fclose(outputFile);
}

int main(int argc, char* argv[]) {
  if (argc == 1) { printUsageAndExit(); }

  if (strcmp(argv[1], "dbshuffle") == 0) {
    if (argc != 4) { printUsageAndExit(); }
    assert(isValidInputFile(argv[2]));
    assert(isValidOutputFile(argv[3]));
    
    shuffleDatabase(argv[2], argv[3]);
  } else {
    printUsageAndExit();
  }
  return 0;
}
