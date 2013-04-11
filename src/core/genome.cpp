#include "core/genome.h"
#include "core/util.h"

using namespace std;

bool readGenome(Genome* g, FILE* inputFilePointer) {
  static char buffer[1010];
  bool iHaveReadSomething = false;

  while(true) {
    fpos_t prevPos;
    fgetpos(inputFilePointer, &prevPos);
    if (!fgets(buffer, sizeof buffer, inputFilePointer)) {
      break;
    }

    if (iHaveReadSomething && buffer[0] == '>') {
      // beginning of a new genome, rollback the file pointer and break
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

bool printGenome(Genome* g, FILE* outputFilePointer, int width) {
  fprintf(outputFilePointer, "%s\n", g->name().c_str());

  for (size_t i = 0; i < g->size(); i += width) {
    for (size_t j = i; j < min(g->size(), i+width); ++j) {
      fprintf(outputFilePointer, "%c", g->data(j));
    }
    fprintf(outputFilePointer, "\n");
  }
}

