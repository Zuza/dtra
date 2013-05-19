#include <cstdio>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "core/Coverage.h"
using namespace std;

const int kNoTests = 20;
const int kN = 30000;
const int randMod = 1000000;

void createIntervals(vector<Interval>* intervals, const int n) {
  for (int i = 0; i < n; ++i) {
    int left = rand() % randMod;
    int right = rand() % randMod;
    if (left > right) { swap(left, right); }

    int value = rand() % randMod;
    intervals->push_back(Interval(left, right, value));
  }
}

bool cmp(const Interval& a, const Interval& b) {
  if (a.value != b.value) return a.value < b.value;
  return a.left < b.left;
}

void drawIntervals(vector<Interval> intervals) {
  sort(intervals.rbegin(), intervals.rend(), cmp);

  for (int row = 0; row < randMod; ++row) {
    vector<char> columns(randMod, '0');
    while (!intervals.empty() && intervals.back().value == row) {
      for (int x = intervals.back().left; x <= intervals.back().right; ++x) {
      	++columns[x];
      }
      intervals.pop_back();
    }
    printf("%3d: ", row);
    for (int col = 0; col < randMod; ++col) {
      if (columns[col] == '0') {
	putchar('.');
      } else {
	putchar(columns[col]);
      }
    } 
    puts("");
  }
}

void printUsageAndExit() {
  puts("test random");
  puts("test <logfile>");
  exit(1);
}

int main(int argc, char* argv[]) {
  if (argc != 2) { printUsageAndExit(); }

  if (strcmp(argv[1], "random") == 0) {
    srand(time(NULL));
    double totalSlowCPU = 0;
    double totalFastCPU = 0;
    
    for (int testNo = 0; testNo < kNoTests; ++testNo) {
      printf("Starting test %d\n", testNo);
      vector<Interval> intervals;
      puts("creating intervals...");
      createIntervals(&intervals, kN);
    
      puts("slow cover...");
      int resultSlow;
      double timerStart = clock();
      coverSlow(NULL, &resultSlow, intervals);
      totalSlowCPU += (clock()-timerStart)/CLOCKS_PER_SEC;
      
      puts("fast cover...");
      int resultFast;
      timerStart = clock();
      cover(NULL, &resultFast, intervals);
      totalFastCPU += (clock()-timerStart)/CLOCKS_PER_SEC;
      
      
      if (resultFast != resultSlow) {
	puts("Test failed! More information available in test.log");
	FILE* log = fopen("test.log", "wt");
	fprintf(log, "Slow result: %d\n", resultSlow);
	fprintf(log, "Fast result: %d\n", resultFast);
	for (size_t i = 0; i < intervals.size(); ++i) {
	  fprintf(log, "[%d %d]: %d\n", intervals[i].left, intervals[i].right,
		  intervals[i].value);
	}
	printf("Slow result: %d\n", resultSlow);
	printf("Fast result: %d\n", resultFast);
	drawIntervals(intervals);
	return 1;
      } else {
	printf("Test %d passed\n", testNo);
      }
    }
    printf("Total time slow: %0.3lf\n", totalSlowCPU);
    printf("Total time fast: %0.3lf\n", totalFastCPU);
  } else {
    FILE* log = fopen(argv[1], "rt");
    char buffer[1000];
    vector<Interval> intervals;
    while (fgets(buffer, sizeof buffer, log)) {
      if (buffer[0] == '[') {
	int left, right, value;
	sscanf(buffer, "[%d %d]: %d", &left, &right, &value);
	intervals.push_back(Interval(left, right, value));
      }
    }
    fclose(log);

    for (size_t i = 0; i < intervals.size(); ++i) {
      printf("[%d %d]: %d\n", intervals[i].left, intervals[i].right,
    	      intervals[i].value);
    }

    puts("slow cover...");
    int resultSlow;
    coverSlow(NULL, &resultSlow, intervals);
    printf("Slow result: %d\n", resultSlow);
    
    puts("fast cover...");
    int resultFast;
    cover(NULL, &resultFast, intervals);
    printf("Fast result: %d\n", resultFast);

    drawIntervals(intervals);
  }
  
  return 0;
}
