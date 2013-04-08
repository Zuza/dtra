#include "lis.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <vector>
using namespace std;

//#define DEBUG_OUTPUT
const int vectorSize = 20000;
const int noTests = 100;
const int maxElement = 1000;

size_t calcLongestIncreasingSubsequenceLength(const vector<size_t>& elements) {
  size_t n = elements.size();
  vector<size_t> dpTable(n);

  for (size_t i = 0; i < n; ++i) {
    dpTable[i] = 1;
    for (size_t j = 0; j < i; ++j) {
      if (elements[j] < elements[i] && dpTable[j]+1 > dpTable[i]) {
	dpTable[i] = dpTable[j]+1;
      }
    }
  }
  
  return *max_element(dpTable.begin(), dpTable.end());
}
				      
void createRandomElements(vector<size_t>* elements, const size_t len) {
  elements->clear();

  for (size_t i = 0; i < len; ++i) {
    size_t elem = rand() % maxElement;
    elements->push_back(elem);
  }
}

double totalSlow = 0;
double totalFast = 0;

void runRandomTest() {
  vector<size_t> elements;
  vector<size_t> result;

  double startRandom = clock();
  createRandomElements(&elements, vectorSize);

  #ifdef DEBUG_OUTPUT
  printf("random: %lf\n", (clock()-startRandom)/CLOCKS_PER_SEC);
  #endif

  double startSlow = clock();
  size_t slowLen = calcLongestIncreasingSubsequenceLength(elements);
  double slowTime = (clock()-startSlow)/CLOCKS_PER_SEC;
  totalSlow += slowTime;

  #ifdef DEBUG_OUTPUT
  printf("slow: %lf\n", slowTime);
  #endif
  
  double startFast = clock();
  calcLongestIncreasingSubsequence(&result, elements);
  double fastTime = (clock()-startFast)/CLOCKS_PER_SEC;
  totalFast += fastTime;

  #ifdef DEBUG_OUTPUT
  printf("fast: %lf\n", fastTime);
  puts("");
  #endif

  assert(result.size() == slowLen);

  // jos provjeri je li dobra sekvenca
  for (size_t i = 1; i < result.size(); ++i) {
    size_t a = result[i-1];
    size_t b = result[i];
    assert(elements[a] < elements[b]);
  }
}

int main(void) {
  srand(time(NULL));

  for (int testId = 0; testId < noTests; ++testId) {
    runRandomTest();
  }

  printf("total slow: %lf\n", totalSlow);
  printf("total fast: %lf\n", totalFast);
  return 0;
}
