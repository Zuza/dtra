#include "BufferedBinaryReader.h"
#include <cstring>

void BufferedBinaryReader::readUnsigned64(unsigned long int* result) {
  assert(sizeof(*result) == 8);
  fread(result, sizeof(*result), 1, inputFile_);
}

void BufferedBinaryReader::readUnsigned32(unsigned int* result) {
  assert(sizeof(*result) == 4);
  fread(result, sizeof(*result), 1, inputFile_);
}

void BufferedBinaryReader::readSigned32(int* result) {
  assert(sizeof(*result) == 4);
  fread(result, sizeof(*result), 1, inputFile_);
}
