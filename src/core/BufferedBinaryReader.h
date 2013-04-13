#ifndef BUFFERED_BINARY_READER
#define BUFFERED_BINARY_READER

#include <cassert>
#include <cstdio>
#include <vector>

class BufferedBinaryReader {
 public:
  BufferedBinaryReader(FILE* inputFile) : 
    inputFile_(inputFile) {
    assert(inputFile_);
  }
  
  void readUnsigned64(unsigned long int* result);
  void readUnsigned32(unsigned int* result);
  void readSigned32(int* result);

  template<typename T> void readVector(std::vector<T>& v) {
    v.clear();
    size_t sz; readUnsigned64(&sz);
    v.resize(sz);
    fread(&v[0], sizeof(T), sz, inputFile_);
  }

 private:
  FILE* inputFile_;
};

#endif
