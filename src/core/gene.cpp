#include "core/gene.h"
#include "core/util.h"

using namespace std;

bool readGene(Gene* g, FILE* inputFilePointer) {
  static char buffer[1010];
  bool iHaveReadSomething = false;
  g->clear();

  size_t nameLen = 0;
  size_t dataLen = 0;


  while(true) {
    fpos_t prevPos;
    fgetpos(inputFilePointer, &prevPos);
    if (!fgets(buffer, sizeof buffer, inputFilePointer)) {
      break;
    }

    if (iHaveReadSomething && buffer[0] == '>') {
      // beginning of a new gene, rollback the file pointer and break
      fsetpos(inputFilePointer, &prevPos);
      break;
    }

    iHaveReadSomething = true;
    string tmp = buffer;
    if (buffer[0] == '>') g->name_ = trim(tmp);
    else g->data_ += trim(tmp);
  }

  return iHaveReadSomething;
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

