/**
 * Database file is processed block by block. Every client gets its own
 * fraction of reads for processing.
 */

#include <unistd.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <iostream>
#include <map>

#include <gflags/gflags.h>

#include "core/ThreadPool.h"
#include "core/database.h"
#include "core/mapping.h"
#include "core/read.h"
#include "core/util.h"

#include <mpi.h>
using namespace std;

DEFINE_int32(kmer_len, 20, "Kmer length");

// readovi se citaju sa stdin-a i salju na stdout
void printUsageAndExit() {
  fprintf(stderr, "stats fasta <fastaPath>\n");
  exit(1);
}

void geneStats(const Gene& g) {
  printf("%s\n", g.name());

  map<string, int> kmer_count;
  for (size_t i = 0; i + FLAGS_kmer_len <= g.dataSize(); ++i) {
    string kmer(g.data()+i, g.data()+i+FLAGS_kmer_len);
    bool ok = true;
    for (int j = 0; j < FLAGS_kmer_len; ++j) {
      kmer[j] = toupper(kmer[j]);
      if (kmer[j] == 'N') {
	ok = false;
      }
    }
    
    if (ok) {
      ++kmer_count[kmer];
    }
  }

  int mini = 1000000000;
  int maxi = -1000000000;
  double mean = 0;
  double sum = 0;
  int n = kmer_count.size();
  for (auto it : kmer_count) {
    mini = min(mini, it.second);
    maxi = max(maxi, it.second);
    mean += it.second;
    sum += it.second;
  }
  mean /= n;
  double stddev = 0;
  for (auto it : kmer_count) {
    stddev += pow(it.second - mean, 2);
    if (it.second == maxi && maxi > mini) {
      printf("%s\n", it.first.c_str());
    }
  }
  stddev = sqrt(stddev/n);

  printf("mini=%d, maxi=%d, len=%d\n", mini, maxi, g.dataSize());
  printf("mean=%lf, stddev=%lf\n", mean, stddev);
  printf("sum=%lf, mean_p=%0.20lf, stddev_p=%0.20lf\n", 
	 sum, mean/sum, stddev/sum);
}

void fastaStats(string fastaInputPath) {
  FILE* finfile = fopen(fastaInputPath.c_str(), "rt");
    
  for (Gene g; readGene(&g, finfile); ) {
    geneStats(g);
  }
  
  fclose(finfile);
}

int main(int argc, char* argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (argc <= 1) {
    printUsageAndExit();
  }
  string command = argv[1];

  if (command == "fasta") {
    if (argc != 3) {
      printUsageAndExit();
    } else {
      assert(isValidInputFile(argv[2]));
    }

    fastaStats(argv[2]);
  } else {
    printUsageAndExit();
  }


  return 0;
}
