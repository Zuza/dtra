// Utility functions and classes for gene assembly.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_UTIL
#define MAPPER_UTIL

#include <ctime>

#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <tr1/unordered_map>

#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

using namespace std;
using namespace std::tr1;

// Containter iterator.
#define ITERATE(x, container) for (__typeof(container.begin()) x = container.begin(); x != container.end(); ++x)

// String splitter function.
inline std::vector<std::string> &Split(
  const std::string &s, 
  char delim, 
  std::vector<std::string> &elems) {

  std::istringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

inline std::vector<std::string> Split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  return Split(s, delim, elems);
}

// Retuns a reverse complement of a given string.
inline std::string ReverseComplement(std::string s) {
  std::string ret;
  ret.reserve(s.size()+1);
  for (int i = (int)s.size()-1; i >= 0; --i) {
    if (s[i] == 'A') ret += "T";
    if (s[i] == 'T') ret += "A";
    if (s[i] == 'C') ret += "G";
    if (s[i] == 'G') ret += "C";
    if (s[i] == 'N') ret += "N";
  }
  return ret;
}

// Simple timer class.
struct Timer {
  Timer(string name) : timer_name_(name){
    start_ = clock();
  }

  void Stop() {
    end_ = clock();
  }

  void Print() {
    printf("%s %lf\n", timer_name_.c_str(), (end_-start_)/CLOCKS_PER_SEC);
  }

  void StopAndPrint() {
    Stop();
    Print();
  }

  string timer_name_;
  double start_;
  double end_;
};

inline int PairSize(const pair<int, int>& p) {
  return p.second - p.first+1;
}

// Class for counting various statistical data.
// Different instances can be mergetd with += operator.
class Counters {
public:
  Counters() {}

  void Add(string key, int val=1) {
    cnt_[key] += val;
  }

  void Print() {
    int total = 0;

    ITERATE(x, cnt_) {
      total += x->second;
    }

    printf("Total: %d\n", total);
    ITERATE(x, cnt_) {
      printf("%s: %d (%lf)\n", x->first.c_str(), x->second, 1.0*x->second/total);
    }
  }

  Counters& operator += (const Counters& a) {
    ITERATE(x, a.cnt_) {
      cnt_[x->first] += x->second;
    }
    return *this;
  }

private:
  map<string, int> cnt_;
};

// Function takes all the candidate positions and merges them together
// to form larger intervals in case any two overlap.
//
// This is to be used with local or semiglobal alignment to eliminate
// excess calculations (on the overlaps).
inline void GroupPositions(const unordered_map<int, int>& candidates,
                    int read_size,
                    vector< pair<int, int> >& res) {

  vector<int> positions;
  positions.reserve(candidates.size());

  // Add all positions and sort them so we can sweep them.
  ITERATE (c, candidates) {
    positions.push_back(c->first);
  }
  sort(positions.begin(), positions.end());
 
  // Last interval in sweeping.
  // Care, hoho! The format is [interval.first, interval_second>.
  pair<int, int> interval(make_pair(0, 0));
 
  ITERATE (p, positions) {
    if ((*p) >= interval.second) {
      if (interval.second) {
        res.push_back(interval);
      }
      interval.first = (*p);
    }
    interval.second = (*p) + read_size;
  }
 
  if (interval.second) {
    res.push_back(interval);
  }
}

// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

#endif  // NAPPER_UTIL
