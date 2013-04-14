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
  printf("client index <database location> <index output directory>\n");
  printf("-- creates the index and dumps its parts in the output dir\n");
  puts("");
  printf("client solve <database location> <index directory> <reads input file>\n");
  printf("-- solves with precomputed index\n");
  exit(1);
}

void inputReads(vector<shared_ptr<Read> >* reads, 
		const string& readsPath, const int limit=-1) {
  int noReads = 0;
  Read tmp;

  FILE* readsIn = fopen(readsPath.c_str(), "rt");
  assert(readsIn);

  for ( ; (limit==-1 || noReads<limit) && tmp.read(readsIn); ++noReads) {
    shared_ptr<Read> newRead(new Read());
    *newRead = tmp;
    reads->push_back(newRead);
  }

  fclose(readsIn);
}


const int seedLen = 16;

void createIndex(string databasePath, string indexFilePath) {
  Database db(databasePath, indexFilePath, seedLen, false);
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
    fprintf(stderr, "Duljine gena [%lu, %lu]\n", minLen, maxLen);
    fprintf(stderr, "U ovom bloku je bilo %d gena.\n", 
            (int)db.getCurrentBlockNoGenes());
    fprintf(stderr, "Prosjecna duljina gena je %lf.\n",
            db.getAverageGeneLength());
    //    checksum = checksum * 10007 + db.checksum();
  }

  //  fprintf(stderr, "DB checksum: %llu\n", checksum);
}

MappingResult solveRead(shared_ptr<Index> idx, shared_ptr<Read> read) {
  MappingResult result;
  performMapping(&result, idx, read);
  return result;
}

void solveReads(const string& databasePath, const string& indexFolderPath, 
		vector<shared_ptr<Read> >& reads) {
  Database db(databasePath, indexFolderPath, seedLen, true);

  for (int indexNo = 0; indexNo < db.getIndexFilesCount(); ++indexNo) {
    shared_ptr<Index> activeIndex = db.readIndexFile(indexNo);

    ThreadPool pool(8); // TODO: ovo staviti ili da automatski detektira
                        // ili da bude jednako broju jezgara na clusteru 
                        // (to je 8)

    vector<future<MappingResult> > results;
    for (int i = 0; i < reads.size(); ++i) {
      results.push_back(pool.enqueue<MappingResult>([i, &activeIndex, &reads] {
	    return solveRead(activeIndex, reads[i]);
	  }));
    } 

    // TODO: skupiti rezultate
  }
}

int main(int argc, char* argv[]) {
  if (argc <= 1) {
    printUsageAndExit();
  }
  string command = argv[1];

  if (command == "index") {
    if (argc != 4) {
      printUsageAndExit();
    }
    createIndex(argv[2], argv[3]);
  } else if (command == "test") {
    // anything here is temporary and can be deleted at any time
    if (argc != 4) {
      printUsageAndExit();
    }
    Database db(argv[2], argv[3], seedLen, true);    
    printf("db count = %d\n", db.getIndexFilesCount());
    shared_ptr<Index> ptr_index = db.readIndexFile(0);
    vector<pair<unsigned int, unsigned int> > ret;
    ptr_index->getPositions(&ret, 1516545110u);
    
    vector<shared_ptr<Gene> >& genes = db.getGenes(); // holds the last loaded index

    for (int i = 0; i <  (int)ret.size(); ++i) {
      printf("(%d, %d)\n", ret[i].first, ret[i].second);
      printf("seq = %s\n", genes[ret[i].first]->data().substr(ret[i].second, 70).c_str());
    }
    
    
  } else if (command == "solve") {
    if (argc != 5) {
      printUsageAndExit();
    }
    vector<shared_ptr<Read> > reads;
    inputReads(&reads, argv[4]);
    solveReads(argv[2], argv[3], reads);
  } else {
    printUsageAndExit();
  }


  return 0;
}
