#include "core/lis.h"
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <mutex>
#include <utility>
#include <vector>
using namespace std;

namespace {

class FenwickMax {
public:
  FenwickMax(int n) {
    elements_ = vector<pair<int, int> > (n+1);
  }
  
  void update(int pos, const pair<int, int>& val) {
    ++pos;

    for ( ; pos < elements_.size(); pos += lobit(pos)) {
      elements_[pos] = max(elements_[pos], val);
    }
  }
  
  pair<int, int> get(int pos) {
    ++pos;

    pair<int, int> ret;
    for ( ; pos > 0; pos -= lobit(pos)) {
      ret = max(ret, elements_[pos]);
    }
    return ret;
  }
  
private:
  int lobit(const int& a) { return a&-a; }
  
private:
  std::vector<std::pair<int, int> > elements_;
};
  
void reconstructLIS(vector<int>* result,
		    int last,
		    const vector<int>& reconstructionTable) {
  result->clear();
  result->reserve(reconstructionTable.size());
  for ( ; last != -1; last = reconstructionTable[last]) {
    result->push_back(last);
  }
  reverse(result->begin(), result->end());
}

};


void calcLongestIncreasingSubsequence(
    vector<int>* result,
    const vector<pair<int, int> >& elements) {
  int n = elements.size();
  if (n == 0) {
    return;
  }

  // for (int i = 1; i < n; ++i) {
  //   assert(elements[i-1].first <= elements[i].first);
  // }

  vector<int> dpTable(n, 0);
  vector<int> reconstructionTable(n, -1);

  int maxSecond = -1;
  for (int i = 0; i < elements.size(); ++i) {
    //assert(elements[i].second >= 0);
    maxSecond = max(maxSecond, elements[i].second);
  }

  // OPASKA: ukoliko je potrebno ova se funkcija moze ubrzati
  // tako da se u 'elements' na ulazi second clan sazme->
  // u tom slucaju potrebno ga je kod rekonstrukcije vratiti
  // na originalnu vrijednost
  FenwickMax fm(maxSecond);

  for (size_t i = 0; i < n; ++i) {
    std::pair<int, int> best = fm.get(elements[i].second-1);
    if (best.first) {
      dpTable[i] = best.first+1;
      reconstructionTable[i] = best.second;
    } else {
      dpTable[i] = 1;
    }
    fm.update(elements[i].second, make_pair(dpTable[i], i));
  }

  int last = max_element(dpTable.begin(), dpTable.end())-dpTable.begin();
  reconstructLIS(result, last, reconstructionTable);
}
