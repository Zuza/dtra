#include "RankedBitmap.hpp"
#include "WaveletTree.hpp"
#include "FmIndex.hpp"

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
  const int data_len = 123123;

  const int num_queries = 1e3;
  const int query_len = 6;

  ////////////////////

  static char data[data_len+123];
  static char orig_data[data_len+123];

  mt19937 mt;
  for (int i = 0; i < data_len; ++i)
    data[i] = alpha[mt() % alpha_sz];
  data[data_len] = 0;
  strncpy(orig_data, data, data_len);
  orig_data[data_len] = 0;

  FmIndex fmindex(data, data_len, alpha, alpha_sz);

  auto slow_count_substrings = [data_len](const char* query, const int query_len, vector<ullint>& results) {
    results.clear();

    ullint ret = 0;
    for (int i = 0; i < data_len; ++i) {
      if (strncmp(orig_data+i, query, query_len) == 0) {
        ++ret;
        results.push_back(i);
      }
    }
    return ret;
  };

  for (int i_q = 0; i_q < num_queries; ++i_q) {
    static char query[query_len + 123];

    for (int i = 0; i < query_len; ++i)
      query[i] = alpha[mt() % alpha_sz];
    query[query_len] = 0;

    static vector<ullint> results, results_brute; results.clear(); results_brute.clear();

    ullint cnt = fmindex.count_substrings(query, query_len);
    ullint cnt_brute = slow_count_substrings(query, query_len, results_brute);

    fmindex.get_substring_pos(results, query, query_len);
    sort(results.begin(), results.end());
    sort(results_brute.begin(), results_brute.end());

    printf("cnt = %llu (brute = %llu)\n", cnt, cnt_brute);
    assert(cnt == cnt_brute);
    assert(results == results_brute);
  }

  return 0; 
}
