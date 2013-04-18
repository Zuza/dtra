#include "core/lis.h"
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <mutex>
#include <utility>
#include <vector>
using namespace std;

namespace {

  // TODO: mozda ce mi intovi tu biti dovoljni umjesto
  // size_t
class FenwickMax {
public:
  FenwickMax(size_t n) {
    elements_ = vector<pair<size_t, size_t> > (n+1);
  }
  
  void update(size_t pos, const pair<size_t, size_t>& val) {
    ++pos;

    for ( ; pos < elements_.size(); pos += lobit(pos)) {
      elements_[pos] = max(elements_[pos], val);
    }
  }
  
  pair<size_t, size_t> get(int pos) {
    ++pos;

    pair<size_t, size_t> ret;
    for ( ; pos > 0; pos -= lobit(pos)) {
      ret = max(ret, elements_[pos]);
    }
    return ret;
  }
  
private:
  size_t lobit(const size_t& a) { return a&-a; }
  
private:
  // napomena: mislim da su intovi dovoljni
  std::vector<std::pair<size_t, size_t> > elements_;
};
  
// TODO: makni rekurziu
void reconstructLIS(vector<int>* result,
		    size_t last,
		    const vector<size_t>& reconstructionTable) {
  if (reconstructionTable[last] != -1) {
    reconstructLIS(result, reconstructionTable[last], reconstructionTable);
  }
  result->push_back(last);
}

};


void calcLongestIncreasingSubsequence(
    vector<int>* result,
    const vector<pair<int, int> >& elements) {
  size_t n = elements.size();

  // TODO: ili maknuti ili ostaviti samo u debug modu
  for (int i = 1; i < n; ++i) {
    assert(elements[i-1].first <= elements[i].first);
  }

  vector<size_t> dpTable(n, 0);
  vector<size_t> reconstructionTable(n, -1);

  int maxSecond = -1;
  for (int i = 0; i < elements.size(); ++i) {
    // TODO: maknuti
    assert(elements[i].second >= 0);

    maxSecond = max(maxSecond, elements[i].second);
  }

  FenwickMax fm(maxSecond);

  for (size_t i = 0; i < n; ++i) {
    std::pair<size_t, size_t> best = fm.get(elements[i].second-1);
    if (best.first) {
      dpTable[i] = best.first+1;
      reconstructionTable[i] = best.second;
    } else {
      dpTable[i] = 1;
    }
    fm.update(elements[i].second, make_pair(dpTable[i], i));
  }

  size_t last = max_element(dpTable.begin(), dpTable.end())-dpTable.begin();
  reconstructLIS(result, last, reconstructionTable);
}
