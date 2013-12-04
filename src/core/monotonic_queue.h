#ifndef MONOTONIC_QUEUE
#define MONOTONIC_QUEUE

#include <cassert>

#include <deque>
#include <functional>
#include <queue>

template<class T>
class MonotonicQueue {
 public:
  void push(const T& val) {
    order.push(val);
    while (!monotonic_vals.empty() && monotonic_vals.back() < val) {
      monotonic_vals.pop_back();
    }
    monotonic_vals.push_back(val);
  }

  void pop() {
    assert(!order.empty());
    assert(!monotonic_vals.empty());

    if (monotonic_vals.front() == order.front()) {
      monotonic_vals.pop_front();
    }
    order.pop();
  }

  const T& front() {
    return order.front();
  }

  const T& max() {
    return monotonic_vals.front();
  }

  int size() {
    return (int)order.size();
  }

  bool empty() {
    return order.empty();
  }

 private:
  std::queue<T> order;
  std::deque<T> monotonic_vals;
};

#endif
