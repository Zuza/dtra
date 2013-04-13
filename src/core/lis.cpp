#include "core/lis.h"
#include <cstdio>
#include <algorithm>
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
  
void reconstructLIS(vector<size_t>* result,
		    size_t last,
		    const vector<size_t>& reconstructionTable) {
  // printf("%d %d\n", last, reconstructionTable[last]);
  if (reconstructionTable[last] != -1) {
    reconstructLIS(result, reconstructionTable[last], reconstructionTable);
  }
  result->push_back(last);
}

};


void calcLongestIncreasingSubsequence(vector<size_t>* result,
				      const vector<size_t>& elements) {
  size_t n = elements.size();
  std::vector<size_t> dpTable(n, 0);
  std::vector<size_t> reconstructionTable(n, -1);

  FenwickMax fm(*max_element(elements.begin(), elements.end()));

  for (size_t i = 0; i < n; ++i) {
    //printf("%d\n", elements[i]);
    std::pair<size_t, size_t> best = fm.get(elements[i]-1);
    if (best.first) {
      dpTable[i] = best.first+1;
      reconstructionTable[i] = best.second;
    } else {
      dpTable[i] = 1;
    }
    fm.update(elements[i], make_pair(dpTable[i], i));
  }

  size_t last = max_element(dpTable.begin(), dpTable.end())-dpTable.begin();
  reconstructLIS(result, last, reconstructionTable);
}
