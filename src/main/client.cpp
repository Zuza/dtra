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
#include <iostream>

#include "core/ThreadPool.h"
#include "core/database.h"
#include "core/read.h"
#include "core/mapping.h"
using namespace std;

// readovi se citaju sa stdin-a i salju na stdout
void printUsageAndExit() {
  printf("client index <database location> <index output location>\n");
  printf("-- creates index and dumps it into a file\n");
  puts("");
  printf("client solve <database location> <index file location>\n");
  printf("-- solves with precomputed index\n");
  exit(1);
}

void inputReads(vector<shared_ptr<Read> >& reads, const int limit) {
  int noReads = 0;
  Read tmp;

  for ( ; noReads < limit && tmp.readFromStdin(); ++noReads) {
    shared_ptr<Read> newRead(new Read());
    *newRead = tmp;
    reads.push_back(newRead);
  }
}

MappingResult solveRead(Database& db, shared_ptr<Read> read) {
  static mutex m;
  MappingResult result;
  performMapping(&result, db, read);
  // m.lock();
  // cout << "Begin thread # " << std::this_thread::get_id() << endl;
  // read->print();
  // cout << "End thread # " << std::this_thread::get_id() << endl;
  // m.unlock();
  return result;
}

void processReads(Database& db, vector<shared_ptr<Read> >& reads) {
  ThreadPool pool(8); // TODO: ovo staviti ili da automatski detektira
                      // ili da bude jednako broju jezgara na clusteru (to je 8)

  vector<future<MappingResult> > results;
  for (int i = 0; i < reads.size(); ++i) {
    results.push_back(pool.enqueue<MappingResult>([i, &db, &reads] {
  	  return solveRead(db, reads[i]);
  	}));
  } 
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    printUsageAndExit();
  }

  string command = argv[1];
  string databasePath = argv[2];
  string indexFilePath = argv[3];

  assert(command == "index" || command == "solve");

  vector<shared_ptr<Read> > reads;
  if (command == "solve") {
    inputReads(reads, 6);
  }

  Database db(databasePath, indexFilePath, 
	      16,
	      command == "index");

  size_t minLen = 1000000000000LL;
  size_t maxLen = 0;
  size_t totalRead = 0;

  unsigned long long checksum = 0;

  for (int blockNumber = 0; db.readNextBlock(); ++blockNumber) {
    if (blockNumber >= 10) break; // TODO: makni limit
    size_t byteLen = db.getCurrentBlockNoBytes();

    if (command == "solve") {
      processReads(db, reads);
    }

    totalRead += byteLen;
    db.getMinMaxGeneLength(&minLen, &maxLen);
    fprintf(stderr, "Read %0.3lf Gb.\n", totalRead/1e9);
    fprintf(stderr, "U ovom bloku je bilo %d gena.\n", 
	    (int)db.getCurrentBlockNoGenes());
    fprintf(stderr, "Prosjecna duljina gena je %lf.\n",
	    db.getAverageGeneLength());
    checksum = checksum * 10007 + db.checksum();
  }
  fprintf(stderr, "Min gene length: %lld bytes.\n", (long long)minLen);
  fprintf(stderr, "Max gene length: %lld bytes.\n", (long long)maxLen);
  fprintf(stderr, "DB checksum: %llu\n", checksum);

  return 0;
}
