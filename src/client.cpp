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

#include "ThreadPool.h"
#include "database.h"
#include "read.h"
#include "mapping.h"
using namespace std;

// readovi se citaju sa stdin-a i salju na stdout
void printUsageAndExit() {
  printf("client index <database location> <index output directory>\n");
  printf("-- creates the index and dumps its parts in the output dir\n");
  puts("");
  printf("client solve <database location> <index directory>\n");
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

// MappingResult solveRead(Database& db, shared_ptr<Read> read) {
//   static mutex m;
//   MappingResult result;
//   //  performMapping(&result, db, read);
//   // m.lock();
//   // cout << "Begin thread # " << std::this_thread::get_id() << endl;
//   // read->print();
//   // cout << "End thread # " << std::this_thread::get_id() << endl;
//   // m.unlock();
//   return result;
// }

// void processReads(Database& db, vector<shared_ptr<Read> >& reads) {
//   ThreadPool pool(8); // TODO: ovo staviti ili da automatski detektira
//                       // ili da bude jednako broju jezgara na clusteru (to je 8)

//   vector<future<MappingResult> > results;
//   for (int i = 0; i < reads.size(); ++i) {
//     results.push_back(pool.enqueue<MappingResult>([i, &db, &reads] {
//   	  return solveRead(db, reads[i]);
//   	}));
//   } 
// }

const int seedLen = 16;

void createIndex(string databasePath, string indexFilePath) {
  Database db(databasePath, indexFilePath, seedLen);
  size_t totalRead = 0;

  //  unsigned long long checksum = 0;

  for (int indexNumber = 0; db.readDbStoreIndex(); ++indexNumber) {
    size_t byteLen = db.getCurrentBlockNoBytes();

    // if (command == "solve") {
    //   processReads(db, reads);
    // }

    totalRead += byteLen;
    size_t minLen, maxLen; db.getMinMaxGeneLength(&minLen, &maxLen);
    fprintf(stderr, "Read %0.3lf Gb.\n", totalRead/1e9);
    fprintf(stderr, "U ovom bloku je bilo %d gena.\n", 
            (int)db.getCurrentBlockNoGenes());
    fprintf(stderr, "Prosjecna duljina gena je %lf.\n",
            db.getAverageGeneLength());
    //    checksum = checksum * 10007 + db.checksum();
  }

  //  fprintf(stderr, "DB checksum: %llu\n", checksum);
}

int main(int argc, char* argv[]) {
  string command = argv[1];

  if (command == "index") {
    if (argc != 4) {
      printUsageAndExit();
    }
    createIndex(argv[2], argv[3]);
  } else if (command == "test") {
    if (argc != 3) {
      printUsageAndExit();
    }
    Database db(argv[2], seedLen);    
    printf("db count = %d\n", db.getIndexFilesCount());
    shared_ptr<Index> ptr_index = db.readIndexFile(0);
    vector<pair<unsigned int, unsigned int> > ret;
    ptr_index->getPositions(&ret, 1516545110u);
    
    for (int i = 0; i <  (int)ret.size(); ++i) {
      printf("(%d, %d)\n", ret[i].first, ret[i].second);
    }
    
  } else if (command == "solve") {
    if (argc != 4) {
      printUsageAndExit();
    }
    assert(false);
    //    vector<shared_ptr<Read> > reads;
    //    inputReads(reads, 6);
  } else {
    printUsageAndExit();
  }


  return 0;
}
