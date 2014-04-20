#include "DnaIndex.hpp"

#include <divsufsort.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>

#include <unordered_map>

using namespace std;

typedef unsigned long long ullint;

void construct_and_save_index() {
  const int num_genes = 3;
  static char gene[][123] = {
    "ACGTAG",
    "AGCAC",
    "GCAATCA",
  };

  DnaIndex dna_index(1e6);

  // for num_genes?
  for (int i = 0; i < num_genes; ++i) {
    dna_index.insert_gene(gene[i], strlen(gene[i]));
  }
  dna_index.create_index();

  FILE* out = fopen("main_dna.index", "w");
  dna_index.serialize(out);
  fclose(out);
}

void read_and_solve() {
  FILE* in = fopen("main_dna.index", "r");
  DnaIndex dna_index(in);
  fclose(in);

  const char query[] = "CA"; const int query_len = strlen(query);
  vector<pair<ullint, ullint> > results;
  dna_index.get_substring_pos(results, query, query_len);

  for (int i = 0; i < (int)results.size(); ++i) {
    printf("gene = %llu position in gene = %llu\n", results[i].first, results[i].second);
  }
}

int main() {
  construct_and_save_index();
  read_and_solve();
  return 0; 
}
