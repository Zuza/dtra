/**
 * Database file is processed block by block. Every client gets its own
 * fraction of reads for processing.
 */

#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <string>
#include "database.h"
using namespace std;

void usage() {
  printf("client index <database location> <index output location>\n");
  printf("-- creates index and dumps it into a file\n");
  puts("");
  printf("client solve <database location> <index file location>\n");
  printf("-- solves with precomputed index\n");
  printf("-- if <index file location> == stdin, index is read from stdin\n");
  exit(1);
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    usage();
  }

  string command = argv[1];
  string databasePath = argv[2];
  string indexFilePath = argv[3];

  assert(command == "index" || command == "solve");

  Database db(databasePath, indexFilePath, 
	      command == "index");

  size_t maxByteLen = 0;
  size_t totalRead = 0;

  for (int blockNumber = 0; db.readNextBlock(); ++blockNumber) {
    //if (blockNumber >= 30) break; // TODO: makni

    size_t byteLen = db.getCurrentBlockNoBytes();

    if (command == "solve") {
      // TODO: ovdje sad radim nesto s readovima
    }
    totalRead += byteLen;
    maxByteLen = max(maxByteLen, byteLen);
    fprintf(stderr, "Read %0.3lf Gb.\n", totalRead/1e9);
    //fprintf(stderr, "Checksum: %llu\n", db.checksum());
  }
  fprintf(stderr, "Max gene length: %lld bytes.\n", (long long)maxByteLen);

  return 0;
}
