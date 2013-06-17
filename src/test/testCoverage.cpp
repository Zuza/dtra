#include <cassert>
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
const int randMod = 3000;

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
  if (kN > 10 || randMod > 20) {
    printf("Example too big, drawing ommited.\n");
    return;
  }

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

void testReconstruction(const vector<Interval>& reconstruction,
			int result,
			const vector<Interval>& intervals,
			const string& debugPrefix) {
  int reconstructionResult = 0;
  for (size_t i = 0; i < reconstruction.size(); ++i) {
    if (i) { assert(reconstruction[i-1].value < 
		    reconstruction[i].value); }
    assert(0 <= reconstruction[i].left);
    assert(reconstruction[i].left <= reconstruction[i].right);
    
    reconstructionResult += 
      reconstruction[i].right -
      reconstruction[i].left + 1;
  }
  
  if (reconstructionResult != result) {
    printf("%s reconstruction not working! :(\n", debugPrefix.c_str());
    drawIntervals(intervals);
    exit(1);
  } else {
    printf("%s reconstruction OK! :)\n", debugPrefix.c_str());
  }
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
      vector<Interval> reconstructionSlow;
      coverSlow(&reconstructionSlow, &resultSlow, intervals);
      totalSlowCPU += (clock()-timerStart)/CLOCKS_PER_SEC;
      
      puts("fast cover...");
      int resultFast;
      timerStart = clock();
      vector<Interval> reconstructionFast;
      cover(&reconstructionFast, &resultFast, intervals);
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
	// TEST SPORE REKONSTRUKCIJE
	testReconstruction(reconstructionSlow, resultSlow, intervals, "slow");

	// TEST BRZE REKONSTRUKCIJE
	//testReconstruction(reconstructionFast, resultFast, intervals, "fast");

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
