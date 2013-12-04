#ifndef FENWICK
#define FENWICK

#include <cassert>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>


template<class T>
class FenwickMax {
public:
  FenwickMax(int n) {
    elements_ = std::vector<T> (n+1, T());
  }
  
  void update(int pos, const T& val) {
    ++pos;

    for ( ; pos < elements_.size(); pos += lobit(pos)) {
      elements_[pos] = std::max(elements_[pos], val);
    }
  }
  

  T get(int pos) {
    ++pos;

    T ret = T();
    for ( ; pos > 0; pos -= lobit(pos)) {
      ret = std::max(ret, elements_[pos]);
    }
    
    //print2(ret);
    return ret;
  }

  // TODO: sluzi za debug, slobodno moze biti obrisano
  void print2(const std::pair<int, int>& a)  {}
  void print2(const int& a) { std::cout << a << std::endl; }
  
private:
  int lobit(const int& a) { return a&-a; }
  
private:
  std::vector<T> elements_;
};


#endif
