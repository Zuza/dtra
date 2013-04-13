#ifndef BUFFERED_BINARY_WRITER
#define BUFFERED_BINARY_WRITER

#include <cassert>
#include <cstdio>
#include <vector>

class BufferedBinaryWriter {
 public:
  BufferedBinaryWriter(FILE* outputFile,
		       const int bufferSize = 1<<20) :
  outputFile_(outputFile),
  bufferSize_(bufferSize) {
    assert(outputFile_);
    assert(bufferSize_ >= 50);
    
    bufferPos_ = 0;
    buffer_ = new char[bufferSize_];
  }

  ~BufferedBinaryWriter() {
    flushBuffer(true);
    delete[] buffer_;
  }

  void writeUnsigned64(const unsigned long int a);
  void writeUnsigned32(const unsigned int a);
  void writeSigned32(const int a);

  template<typename T> void writeVector(const std::vector<T>& v) {
    writeUnsigned64(v.size());
    flushBuffer(true);
    fwrite((void*)&v[0], sizeof(T), v.size(), outputFile_);
  }

  void flushBuffer(bool forceFlush = false);

 private:
  FILE* outputFile_;
  int bufferSize_;
  size_t bufferPos_;
  char* buffer_;
};

#endif
