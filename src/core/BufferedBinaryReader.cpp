#include "core/BufferedBinaryReader.h"
#include <cstring>

void BufferedBinaryReader::updateBuffer() {
  if (bufferPos_+10 >= bufferSize_) {
    int diff = bufferSize_-bufferPos_;
    if (diff>0) {
      memcpy((void*)buffer_, (void*)(buffer_+bufferPos_), diff);
    }
    fread((void*)(buffer_+diff), 1, bufferSize_-diff, inputFile_);
    bufferPos_ = 0;
  }
}

void BufferedBinaryReader::readUnsigned64(unsigned long long* result) {
  *result = *((unsigned long long*)(buffer_+bufferPos_));
  bufferPos_ += 8;
  updateBuffer();
}

void BufferedBinaryReader::readUnsigned32(unsigned int* result) {
  *result = *((unsigned int*)(buffer_+bufferPos_));
  bufferPos_ += 4;
  updateBuffer();
}

void BufferedBinaryReader::readSigned64(long long* result) {
  *result = *((long long*)(buffer_+bufferPos_));
  bufferPos_ += 8;
  updateBuffer();
}

void BufferedBinaryReader::readSigned32(int* result) {
  *result = *((int*)(buffer_+bufferPos_));
  bufferPos_ += 4;
  updateBuffer();
}

void BufferedBinaryReader::readChar(char* result) {
  *result = *(buffer_+bufferPos_);
  bufferPos_ += 1;
  updateBuffer();
}



