#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <string>
#include "core/genome.h"
using namespace std;

inline bool throwCoin(double p) {
  const int precision = 1000;
  return rand() % precision < precision*p;
}

void printUsageAndExit() {
  // TODO: opcionalizirati kroz komandnu liniju
  printf("Usage:\n");
  printf("reducer nt <nt input file>\n");
  printf("reducer wgsim <nt input file>\n");
  exit(1);
}

void reduceNtDatabase(char* ntFilePath) {
  FILE* ntInputFile = fopen(ntFilePath, "rt");

  // TODO: opcionalizirati ovu funkciju
  int koji = 0;
  for (Genome g; readGenome(&g, ntInputFile); ++koji) {
    if (koji % 1000 == 0) {
      fprintf(stderr, "Processed %d genes.\n", koji);
    }
    if (koji >= 700) {
      break;
    }
    //if (throwCoin(0.1)) {
      printGenome(&g);
      //}
  }
  
  fclose(ntInputFile);
}

void createWgsimReads(char* ntFilePath) {
  FILE* ntInputFile = fopen(ntFilePath, "rt");
  // TODO: dati mogucnost da se wgsim opcije postavljaju kroz ovaj program

  for (Genome g; readGenome(&g, ntInputFile); ) {
    FILE* tmpGeneFile = fopen("reducer.gene.tmp", "wt");
    printGenome(&g, tmpGeneFile);

    const string wgsim = "./wgsim ";
    string command;
    assert(system(""));

    fclose(tmpGeneFile);
  }

  fclose(ntInputFile);
}

int main(int argc, char* argv[]) {
  srand(time(NULL));

  if (argc != 3 || strcmp(argv[1], "nt")) {
    printUsageAndExit();
  }

  if (strcmp(argv[1], "nt") == 0) {
    char* filepath = argv[2];
    reduceNtDatabase(filepath);
  } else if (strcmp(argv[1], "wgsim")) {
    char* filepath = argv[2];
    createWgsimReads(filepath);
  }
  return 0;
}
