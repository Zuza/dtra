#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <string>
#include "core/gene.h"
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
  for (Gene g; readGene(&g, ntInputFile); ++koji) {
    if (koji % 1000 == 0) {
      fprintf(stderr, "Processed %d genes.\n", koji);
    }
    if (koji >= 10) {
      break;
    }
    //if (throwCoin(0.1)) {
      printGene(&g);
      //}
  }
  
  fclose(ntInputFile);
}

void createWgsimReads(char* ntFilePath) {
  FILE* ntInputFile = fopen(ntFilePath, "rt");
  assert(ntInputFile);
  
  const string outputReadsBig = "reads.fq";
  system(("rm " + outputReadsBig).c_str());

  for (Gene g; readGene(&g, ntInputFile); ) {
    const string tmpFasta = "reducer.gene.tmp.fa";
    FILE* tmpGeneFile = fopen(tmpFasta.c_str(), "wt");
    printGene(&g, tmpGeneFile);
    fclose(tmpGeneFile);

    const string wgsim = "./wgsim ";
    const int readsPerGene = 15;
    
    const int readLength = min((int)g.size(), 700);
    const string readsOutput1 = "reads.output.tmp.1";
    const string readsOutput2 = "reads.output.tmp.2";
    
    ostringstream command;
    command << wgsim;
    command << " -N " << readsPerGene;
    command << " -1 " << readLength;
    command << " -2 " << readLength;
    command << " " << tmpFasta;
    command << " " << readsOutput1;
    command << " " << readsOutput2;
    command << " > /dev/null";
    system(command.str().c_str());

    system(("cat " + readsOutput1 + " >> " + outputReadsBig).c_str());
    system(("rm " + readsOutput1).c_str());
    system(("rm " + readsOutput2).c_str());
    system(("rm " + tmpFasta).c_str());
  }

  fclose(ntInputFile);
}

int main(int argc, char* argv[]) {
  srand(time(NULL));

  if (argc != 3) {
    printUsageAndExit();
  }

  if (strcmp(argv[1], "nt") == 0) {
    char* filepath = argv[2];
    reduceNtDatabase(filepath);
  } else if (strcmp(argv[1], "wgsim") == 0) {
    char* filepath = argv[2];
    createWgsimReads(filepath);
  } else {
    printUsageAndExit();
  }
  return 0;
}
