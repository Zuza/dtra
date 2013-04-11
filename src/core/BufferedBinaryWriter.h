#ifndef BUFFERED_BINARY_WRITER
#define BUFFERED_BINARY_WRITER

#include <cassert>
#include <cstdio>

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

  void writeUnsigned32(const unsigned int& a);
  void writeSigned32(const int& a);
  void flushBuffer(bool forceFlush = false);

 private:
  FILE* outputFile_;
  int bufferSize_;
  size_t bufferPos_;
  char* buffer_;
};

#endif
