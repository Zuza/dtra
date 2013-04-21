/**
 * Database file is processed block by block. Every client gets its own
 * fraction of reads for processing.
 */

#include <unistd.h>

#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <iostream>
#include <map>

#include <gflags/gflags.h>

#include "core/ThreadPool.h"
#include "core/database.h"
#include "core/index.h"
#include "core/mapping.h"
#include "core/read.h"
using namespace std;

DEFINE_int32(seed_len, 20, "Seed length that is stored/read from the index");
DEFINE_bool(validate_wgsim, false, "Used for simulated tests, if true some statistics is printed on stdout.");

// readovi se citaju sa stdin-a i salju na stdout
void printUsageAndExit() {
  printf("client index <database location> <index output directory>\n");
  printf("client solve <database location> <index directory> <reads input file> <result output file>\n");
  exit(1);
}

void createIndex(string databasePath, string indexFilePath) {
  Database db(databasePath, indexFilePath, FLAGS_seed_len, false);
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
    fprintf(stderr, "Begin processing block %d/%d...\n", 
	    indexNo+1, indexFileCount);
    shared_ptr<Index> activeIndex = db.readIndexFile(indexNo);
    vector<shared_ptr<Gene> >& currentGenes = db.getGenes();

    //http://stackoverflow.com/questions/150355/
    //programmatically-find-the-number-of-cores-on-a-machine
    int numCores = sysconf(_SC_NPROCESSORS_ONLN);
    assert(numCores > 1 && numCores < 100); // sanity check
    ThreadPool pool(numCores-1); // one core for this thread
    
    vector<future<int> > results;
    for (int i = 0; i < reads.size(); ++i) {
      results.push_back(
	 pool.enqueue<int>([i, &currentGenes, &activeIndex, &reads] {
	    return solveRead(currentGenes, activeIndex, reads[i]);
	  }));
    } 

    for (int i = 0; i < reads.size(); ++i) {
      results[i].wait();
    }
  }
}

void printWgsimStatistics(const vector<shared_ptr<Read> >& reads) {
  map<int, int> stats;
  for (int i = 0; i < reads.size(); ++i) {
    shared_ptr<Read> read = reads[i];

    int mappingQuality = read->validateWgsimMapping();
    ++stats[mappingQuality];

    if (mappingQuality == -1) {      
      printf("READ #%04d:\n", i);
      printf("id: %s\n", read->id().c_str());
      printf("data: %s\n", read->data().c_str());
      printf("mappings (%d):\n", (int)read->topMappings().size());
      
      for (int x = 0; x < read->topMappings().size(); ++x) {
        OneMapping mapping = read->topMapping(x);
        printf("score=%lf geneDes=%s genePos=%d isRC=%d\n",
               mapping.score, mapping.geneDescriptor.c_str(), 
	       mapping.genePos, mapping.isRC);
      }
      puts("END READ");
      puts("");
    }
  }

  printf("Total reads: %d\n", (int)reads.size());
  printf("Number of reads not mapped: %d\n", stats[-1]);
  for (map<int, int>::iterator it = stats.begin(); it != stats.end(); ++it) {
    if (it->first != -1) {
      printf("hitova na %d-tom mjestu: %d\n", it->first+1, it->second);
    }
  }
}

void printReads(const vector<shared_ptr<Read> >& reads,
		const string& resultFilePath) {
  FILE* resultOut = fopen(resultFilePath.c_str(), "wt"); 
  // format outputa: read_id,top_aln_num;nucl_id,score,start,stop,strand;...
  assert(resultOut);

  for (int i = 0; i < reads.size(); ++i) {
    shared_ptr<Read> read = reads[i];

    fprintf(resultOut, "%s,%d", read->id().c_str(), 
	    (int)read->topMappings().size());
    
    for (int x = 0; x < read->topMappings().size(); ++x) {
      const OneMapping& onemap = read->topMapping(x);
      
      fprintf(resultOut, ";%s,%lf,%d,%d,%d", onemap.geneDescriptor.c_str(),
      	      onemap.score, onemap.genePos, 
      	      onemap.genePos+read->size(), onemap.isRC);
    }
    fprintf(resultOut, "\n");
  }

  fclose(resultOut);

  if (FLAGS_validate_wgsim) {
    printWgsimStatistics(reads);
  }
}

int main(int argc, char* argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);

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

    Database db(argv[2], argv[3], FLAGS_seed_len, true);    
    printf("db count = %d\n", db.getIndexFilesCount());
    shared_ptr<Index> ptr_index = db.readIndexFile(0);
    vector<pair<unsigned int, unsigned int> > ret;
    ptr_index->getPositions(&ret, 1516545110u);
    
    vector<shared_ptr<Gene> >& genes = db.getGenes(); // holds the last loaded index

    for (int i = 0; i <  (int)ret.size(); ++i) {
      printf("(%d, %d)\n", ret[i].first, ret[i].second);
      printf("seq = %s\n", genes[ret[i].first]->data());
    }
    
  } else if (command == "solve") {
    if (argc != 6) {
      printUsageAndExit();
    }
    Database db(argv[2], argv[3], FLAGS_seed_len, true);    
    vector<shared_ptr<Read> > reads;
    inputReads(&reads, argv[4]); 
    solveReads(db, reads); 
    printReads(reads, argv[5]);
  } else {
    printUsageAndExit();
  }


  return 0;
}
