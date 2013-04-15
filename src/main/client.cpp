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
  printf("client solve <database location> <index directory> <reads input file>\n");
  exit(1);
}

const int kSeedLen = 16;

void createIndex(string databasePath, string indexFilePath) {
  Database db(databasePath, indexFilePath, kSeedLen, false);
  size_t totalRead = 0;

  for (int indexNumber = 0; db.readDbStoreIndex(); ++indexNumber) {
    size_t byteLen = db.getCurrentBlockNoBytes();

    totalRead += byteLen;
    size_t minLen, maxLen; db.getMinMaxGeneLength(&minLen, &maxLen);
    fprintf(stderr, "Read %0.3lf Gb.\n", totalRead/1e9);
    fprintf(stderr, "Duljine gena [%lu, %lu]\n", minLen, maxLen);
    fprintf(stderr, "U ovom bloku je bilo %d gena.\n", 
            (int)db.getCurrentBlockNoGenes());
    fprintf(stderr, "Prosjecna duljina gena je %lf.\n",
            db.getAverageGeneLength());
  }
}

void inputReads(vector<shared_ptr<Read> >* reads, 
		const string& readsPath, const int limit=-1) {
  FILE* readsIn = fopen(readsPath.c_str(), "rt");
  assert(readsIn);

  Read tmp;
  int noReads = 0;
  for ( ; (limit==-1 || noReads<limit) && tmp.read(readsIn); ++noReads) {
    shared_ptr<Read> newRead(new Read());
    *newRead = tmp;
    reads->push_back(newRead);
  }

  fclose(readsIn);
}

int solveRead(vector<shared_ptr<Gene> >& genes, 
	      shared_ptr<Index> idx, shared_ptr<Read> read) {
  performMapping(genes, idx, read);
  return 0;
}

void solveReads(Database& db, 
		vector<shared_ptr<Read> >& reads) {
  int indexFileCount = db.getIndexFilesCount();
  for (int indexNo = 0; indexNo < indexFileCount; ++indexNo) {
    shared_ptr<Index> activeIndex = db.readIndexFile(indexNo);
    vector<shared_ptr<Gene> >& genes = db.getGenes();

    ThreadPool pool(8); // TODO: ovo staviti ili da automatski detektira
                        // ili da bude jednako broju jezgara na clusteru 
                        // (to je 8)

    vector<future<int> > results;
    for (int i = 0; i < reads.size(); ++i) {
      results.push_back(pool.enqueue<int>([i, &genes, &activeIndex, &reads] {
	    return solveRead(genes, activeIndex, reads[i]);
	  }));
    } 

    for (int i = 0; i < reads.size(); ++i) {
      results[i].wait();
    }
  }
}

void printReads(Database& db, const vector<shared_ptr<Read> >& reads) {
  vector<shared_ptr<Gene> >& genes = db.getGenes(); // holds the last loaded index  

  int readsNoMappings = 0;

  for (int i = 0; i < reads.size(); ++i) {
    printf("READ #%04d:\n", i);
    shared_ptr<Read> read = reads[i];
    printf("id: %s\n", read->id().c_str());
    printf("data: %s\n", read->data().c_str());
    printf("mappings (%d):\n", (int)read->topMappings().size());
    if (read->topMappings().size() == 0) {
      ++readsNoMappings;
    }

    for (int x = 0; x < read->topMappings().size(); ++x) {
      OneMapping mapping = read->topMapping(x);
      printf("score=%lf geneId=%d genePos=%d isRC=%d\n",
	     mapping.score, mapping.geneId, mapping.genePos, mapping.isRC);
      // NE VALJA SLJEDECA LINIJA KAD IMAM VISE OD 1 BLOKA
      printf("gene: %s\n", genes[mapping.geneId]->name().c_str());
      printf("segment: %s\n", genes[mapping.geneId]->data(mapping.genePos,
					  mapping.genePos+read->size()).c_str());
    }
    
    puts("END READ");
    puts("");
  }

  printf("Total reads: %d\n", (int)reads.size());
  printf("Number of reads not mapped: %d\n", readsNoMappings);

  for (int i = 0; i < (int)genes.size(); ++i) {
    printGene(genes[i].get());
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

    Database db(argv[2], argv[3], kSeedLen, true);    
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
    Database db(argv[2], argv[3], kSeedLen, true);    
    vector<shared_ptr<Read> > reads;
    inputReads(&reads, argv[4]);
    solveReads(db, reads);
    printReads(db, reads);
  } else {
    printUsageAndExit();
  }


  return 0;
}
