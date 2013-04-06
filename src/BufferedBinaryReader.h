#ifndef BUFFERED_BINARY_READER
#define BUFFERED_BINARY_READER

#include <cassert>
#include <cstdio>

class BufferedBinaryReader {
 public:
  BufferedBinaryReader(FILE* inputFile,
		       const int bufferSize = 1<<20) : 
  inputFile_(inputFile),
  bufferSize_(bufferSize) {
    assert(inputFile_);
    assert(bufferSize_ >= 50);

    bufferPos_ = 0;
    buffer_ = new char[bufferSize_];
    fread((void*)buffer_, 1, bufferSize_, inputFile_);
  }
  
  ~BufferedBinaryReader() {
    delete[] buffer_;
  }

  void readUnsigned64(unsigned long long* result);
  void readUnsigned32(unsigned int* result);
  void readSigned64(long long* result);
  void readSigned32(int* result);
  void readChar(char* result);

 private:
  void updateBuffer();

 private:
  FILE* inputFile_;
  int bufferSize_;
  size_t bufferPos_;
  char* buffer_;
};

#endif
