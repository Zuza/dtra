#include "RankedBitmap.hpp"
#include "WaveletTree.hpp"
#include "FmIndex.hpp"

#include "RankedBitmap.hpp"
#include "WaveletTree.hpp"
#include "FmIndex.hpp"

#include "util.hpp"

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

int main() {
  const char alpha[] = "ACGT"; int alpha_sz = strlen(alpha);
  const int data_len = 123123123;

  const int num_queries = 1e6;
  const int query_len = 20;

  ////////////////////

  int start_time, end_time;

  static char data[data_len+123];

  mt19937 mt;
  for (int i = 0; i < data_len; ++i)
    data[i] = alpha[mt() % alpha_sz];
  data[data_len] = 0;

  start_time = clock();
  FmIndex fmindex(data, data_len, alpha, alpha_sz);
  end_time = clock();

  fprintf(stderr, "index construction time = %.2lf\n", -(start_time - end_time) / double(CLOCKS_PER_SEC));

  ullint tot_results = 0;
  mt.seed(12347);
  start_time = clock();  
  for (int i_q = 0; i_q < num_queries; ++i_q) {
    static char query[query_len + 123];

    for (int i = 0; i < query_len; ++i)
      query[i] = alpha[mt() % alpha_sz];
    query[query_len] = 0;

    static vector<ullint> results;

    fmindex.get_substring_pos(results, query, query_len);
    tot_results += results.size();
  }
  end_time = clock();

  fprintf(stderr, "query time = %.2lf\n", -(start_time - end_time) / double(CLOCKS_PER_SEC));
  fprintf(stderr, "total number of results returned = %llu\n", tot_results);


  ////////////// serialize ////////////////////

  FILE* out = fopen("fmindex.out", "w");
  fmindex.serialize(out);
  fclose(out);

  return 0; 
}
