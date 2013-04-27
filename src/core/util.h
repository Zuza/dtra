// Utility functions and classes for gene assembly.
//
// Authors: Filip Pavetic, Goran Zuzic

#ifndef MAPPER_UTIL
#define MAPPER_UTIL

#include <cctype>
#include <ctime>

#include <algorithm>
#include <functional> 
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <sys/stat.h>
#include <gflags/gflags.h>

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

static inline bool isValidInputFile(const std::string& filePath) {
  FILE* f = fopen(filePath.c_str(), "r");
  bool ok = (f != NULL);
  fclose(f);
  return ok;
}

static inline bool isValidOutputFile(const std::string& filePath) {
  FILE* f = fopen(filePath.c_str(), "w");
  bool ok = (f != NULL);
  fclose(f);
  return ok;
}

static inline unsigned long long getFileSize(const std::string& filePath) {
  struct stat st;
  stat(filePath.c_str(), &st);
  return st.st_size;
}

static inline std::string cstrToString(const char* s, int len) {
  return std::string(s, s + len);
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

#endif  // MAPPER_UTIL
