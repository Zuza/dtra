#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>

#include "core/genome.h"

inline bool throwCoin(double p) {
  const int precision = 1000;
  return rand() % precision < precision*p;
}

void printUsageAndExit() {
  printf("Usage:\n");
  printf("reducer nt <nt input file>\n");
  exit(1);
}

void reduceNtDatabase(char* ntFilePath) {
  FILE* ntInputFile = fopen(ntFilePath, "rt");

  for (Genome g; readGenome(&g, ntInputFile); ) {
    if (throwCoin(1)) {
      printGenome(&g);
    }
  }
  
  fclose(ntInputFile);
}

int main(int argc, char* argv[]) {
  srand(time(NULL));

  if (argc != 3 || strcmp(argv[1], "nt")) {
    printUsageAndExit();
  }

  if (strcmp(argv[1], "nt") == 0) {
    reduceNtDatabase(argv[2]);
  }
  return 0;
}
