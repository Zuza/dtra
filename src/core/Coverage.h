#ifndef COVERAGE
#define COVERAGE

#include <vector>

struct Interval {
  int left, right, value;
  Interval(int left, int right, int value) :
    left(left), right(right), value(value) {}
};

void coverFast(std::vector<Interval>* resultCoverage,
	       int* result,
	       const std::vector<Interval>& intervals);

void coverSlow(std::vector<Interval>* resultCoverage,
	       int* result,
	       std::vector<Interval> intervals);

void cover(std::vector<Interval>* resultCoverage,
	   int* result,
	   const std::vector<Interval>& intervals);

#endif
