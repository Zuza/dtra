#include "BufferedBinaryWriter.h"
#include <cstring>

using namespace std;

void BufferedBinaryWriter::flushBuffer(bool forceFlush) {
  if (forceFlush || bufferPos_+10 >= bufferSize_) {
    fwrite((void*)buffer_, 1, bufferPos_, outputFile_);
    bufferPos_ = 0;
  }
}

void BufferedBinaryWriter::writeUnsigned64(const unsigned long int a) {
  //*((unsigned int*)(buffer_+bufferPos_)) = a;
  assert(sizeof(a) == 8);
  memcpy((void*)(buffer_+bufferPos_), (void*)&a, sizeof(a));
  bufferPos_ += sizeof(a);
  flushBuffer();
}

void BufferedBinaryWriter::writeUnsigned32(const unsigned int a) {
  //*((unsigned int*)(buffer_+bufferPos_)) = a;
  assert(sizeof(a) == 4);
  memcpy((void*)(buffer_+bufferPos_), (void*)&a, sizeof(a));
  bufferPos_ += sizeof(a);
  flushBuffer();
}

void BufferedBinaryWriter::writeSigned32(const int a) {
  //*((unsigned int*)(buffer_+bufferPos_)) = a;
  assert(sizeof(a) == 4);
  memcpy((void*)(buffer_+bufferPos_), (void*)&a, sizeof(a));
  bufferPos_ += sizeof(a);
  flushBuffer();
}
