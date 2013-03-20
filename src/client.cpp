/**
 * Database file is processed block by block. Every client gets its own
 * fraction of reads for processing.
 */

#include <cstdio>
#include <cstring>
#include "database.h"

int main(int argc, char* argv[]) {
  Database db(argv[1]);

  for (size_t totalRead = 0; db.readNextBlock(); ) {
    totalRead += db.getCurrentBlockNoBytes();
    fprintf(stderr, "Read %lld bytes.\n", (long long)totalRead);
  }
  //db.printNames();

  return 0;
}
