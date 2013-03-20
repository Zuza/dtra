#include "genome.h"
#include "util.h"

using namespace std;

bool readGenome(Genome* g, FILE* inputFilePointer) {
  static char buffer[1010];
  bool iHaveReadSomething = false;

  while(true) {
    fpos_t prevPos;
    //printf("%p\n", inputFilePointer);
    fgetpos(inputFilePointer, &prevPos);
    //puts("tu");
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


