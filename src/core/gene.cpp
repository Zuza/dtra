#include "core/gene.h"
#include "core/util.h"
#include "core/bioinf_util.h"

#include <cassert>
#include <cctype>

using namespace std;

bool readGene(Gene* g, FILE* inputFilePointer) {
  static char buffer[1010];
  g->clear();

  size_t nameLen = 0;
  size_t dataLen = 0;
  size_t numGetsCallsName = 0;
  size_t numGetsCallsData = 0;

  fpos_t prevPos;
  fgetpos(inputFilePointer, &prevPos);

  bool doneReadingName = false;
  while (!doneReadingName) {  // read name()
    if (!fgets(buffer, sizeof buffer, inputFilePointer)) {
      assert(dataLen == 0 && nameLen == 0);
      return false;
    } else if (nameLen == 0) {
      assert(buffer[0] == '>');
    }

    int len = strlen(buffer);
    if (buffer[len-1] == '\n') {
      doneReadingName = true;
    }
    while (len > 0 && isspace(buffer[len-1])) {
      buffer[len-1] = 0;
      --len;
    }

    nameLen += len;
    ++numGetsCallsName;
  }
  
  bool doneReadingData = false;
  while (!doneReadingData) {  // read data()
    if (!fgets(buffer, sizeof buffer, inputFilePointer)) {
      break;
    } else if (buffer[0] == '>') {
      break;
    }

    int len = strlen(buffer);
    while (len > 0 && isspace(buffer[len-1])) {
      buffer[len-1] = 0;
      --len;
    }

    dataLen += len;
    ++numGetsCallsData;
  }

  g->name_ = (char*)malloc(nameLen + 1);
  assert(g->name_);
  g->name_len_ = nameLen;

  g->data_ = (char*)malloc(dataLen + 1);
  assert(g->data_);
  g->data_len_ = dataLen;

  size_t ptrName = 0;

  fsetpos(inputFilePointer, &prevPos);
  for (size_t i = 0; i < numGetsCallsName; ++i) { 
    fgets(buffer, sizeof buffer, inputFilePointer);

    int len = strlen(buffer);
    while (len > 0 && isspace(buffer[len-1])) {
      buffer[len-1] = 0;
      --len;
    }

    memcpy(g->name_ + ptrName, buffer, len);
    ptrName += len;
  }

  size_t ptrData = 0;
  for (size_t i = 0; i < numGetsCallsData; ++i) { 
    fgets(buffer, sizeof buffer, inputFilePointer);

    int len = strlen(buffer);
    while (len > 0 && isspace(buffer[len-1])) {
      buffer[len-1] = 0;
      --len;
    }

    memcpy(g->data_ + ptrData, buffer, len);
    ptrData += len;
  }

  g->data_[dataLen] = 0;
  g->name_[nameLen] = 0;

  assert(ptrData == dataLen);
  assert(ptrName == nameLen);

  for (size_t i = 0; i < dataLen; ++i) {
    assert(g->data_[i] >= 'A' && g->data_[i] <= 'Z');
  }

  return dataLen > 0;
}

size_t printGene(Gene* g, FILE* outputFilePointer, int width) {
  size_t printed = 0;
  printed += fprintf(outputFilePointer, "%s\n", g->name());

  for (size_t i = 0; i < g->size(); i += width) {
    for (size_t j = i; j < min(g->size(), i+width); ++j) {
      printed += fprintf(outputFilePointer, "%c", g->data(j));
    }
    printed += fprintf(outputFilePointer, "\n");
  }

  return printed;
}

