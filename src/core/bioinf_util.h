#ifndef BIOINF_UTIL
#define BIOINF_UTIL

#include <cassert>
#include <ctime>
#include <cstdlib>

inline bool isBase(char isit) {
  isit = toupper(isit);

  return isit == 'A' || isit == 'T' || isit == 'G' || isit == 'C';
}

inline int baseToInt(char base, const int nValue = -1) {
  base = toupper(base);

  if (base == 'A') return 0;
  if (base == 'C') return 1;
  if (base == 'T') return 2;
  if (base == 'G') return 3;
  
  // by default, N will return random number [0..3]
  if (nValue == -1) return rand()%4;
  return nValue;
}

inline char intToBase(const int num) {
  if (num == 0) return 'A';
  if (num == 1) return 'C';
  if (num == 2) return 'T';
  if (num == 3) return 'G';
  
  // hopefully this doesn't happen
  assert(0 <= num && num < 4);
  return 0;
}

inline char getBaseComplement(char base) {
  base = toupper(base);

  if (base == 'A') return 'T';
  if (base == 'T') return 'A';
  if (base == 'C') return 'G';
  if (base == 'G') return 'C';
  assert(base == 'N');
  return 'N';
}

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

// translate bases (ie. 'Y' -> C or T), ref: http://www.bioinformatics.org/sms/iupac.html

inline char randBaseToBase(const char base) {
  static int len[26] = {1, 3, 1, 3, 0, 0, 1, 3, 0, 0, 2, 0, 2, 4, 0, 0, 0, 2, 2, 1, 1, 3, 2, 0, 2, 0};
  static char trans[26][5] = {"A",
                              "CGT",
                              "C",
                              "AGT",
                              "",
                              "",
                              "G",
                              "ACT",
                              "",
                              "",
                              "GT",
                              "",
                              "AC",
                              "AGTC",
                              "",
                              "",
                              "",
                              "AG",
                              "GC",
                              "T",
                              "U",
                              "ACG",
                              "AT",
                              "",
                              "CT",
                              ""};

  assert(base >= 'A' && base <= 'Z');

  if (len[base-'A'] == 1) {
    return base;
  }
  
  int t = len[base-'A'];
  assert(t != 0);
  t = rand() % t;
  return trans[base-'A'][t];
}

#endif
