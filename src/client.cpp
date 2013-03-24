/**
 * Database file is processed block by block. Every client gets its own
 * fraction of reads for processing.
 */

#include <cstdio>
#include <cstring>
#include <algorithm>
#include "database.h"
using namespace std;

int main(int argc, char* argv[]) {
  Database db(argv[1]);


  size_t maxByteLen = 0;
  for (size_t totalRead = 0; db.readNextBlock(); ) {
    size_t byteLen = db.getCurrentBlockNoBytes();
    totalRead += byteLen;
    maxByteLen = max(maxByteLen, byteLen);
    fprintf(stderr, "Read %0.3lf Gb.\n", totalRead/1e9);
  }
  fprintf(stderr, "Max gene length: %lld bytes.\n", (long long)maxByteLen);

  return 0;
}
