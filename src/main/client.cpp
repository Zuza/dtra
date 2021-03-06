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

#include "FmIndexWavelet/DnaIndex.hpp"
#include "core/ThreadPool.h"
#include "core/database.h"
#include "core/mapping.h"
#include "core/read.h"
#include "core/util.h"

#include <mpi.h>
using namespace std;

DEFINE_int32(seed_len, 20, "Seed length that is stored/read from the index");
DEFINE_int32(solver_threads, sysconf(_SC_NPROCESSORS_ONLN),
             "Number of threads used by the solver");
DEFINE_string(validate_wgsim, "", "report/full");
DEFINE_int32(no_reads, -1, "Number of reads to process.");
DEFINE_double(confidence, 0.5, "Confidence threshold.h");

DEFINE_bool(calc_edit_distance, true, "Calculate actual edit-distance for read placements?");

// readovi se citaju sa stdin-a i salju na stdout
void printUsageAndExit() {
  printf("client index <database location> <index output directory>\n");
  printf("client solve <database location> <index directory> <reads input file> <result output file>\n");
  printf("mpirun <mpi arguments> client cluster <database location> <index directory> <reads input file> <result output file>\n");
  exit(1);
}

void createIndex(string databasePath, string indexFilePath) {
  Database db(databasePath, indexFilePath, FLAGS_seed_len, false);
  db.readDbStoreIndex(); // creates and stores the index
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
              shared_ptr<DnaIndex> idx,
              const int seedLen,
              shared_ptr<Read> read,
              bool fillEditDistance) {
  performMapping(genes, idx, seedLen, read, fillEditDistance);
  return 0;
}

void solveReads(Database& db, 
                const int seedLen,
                vector<shared_ptr<Read> >& reads,
                bool fillEditDistance) {
  int index_file_count = db.getIndexFilesCount();

  clock_t starting_time;
  TRACE(index_file_count);
  for (int index_no = 0; index_no < index_file_count; ++index_no) {
    starting_time = clock();
    fprintf(stderr, "Processing block %d/%d... ", index_no+1, index_file_count);
    shared_ptr<DnaIndex> activeIndex = db.readIndexFile(index_no);
    vector<shared_ptr<Gene> >& currentGenes = db.getGenes();

    // Find the number of cores on this machine
    // http://stackoverflow.com/questions/150355/programmatically-find-the-number-of-cores-on-a-machine
    int threads = FLAGS_solver_threads;
    assert(threads >= 1 && threads < 100); // sanity check
    ThreadPool pool(threads); // one core for this thread

    vector<future<int> > results;
    for (int i = 0; i < reads.size(); ++i) {
      results.push_back(
        pool.enqueue<int>([i, seedLen, &currentGenes, &activeIndex, &reads, fillEditDistance] {
            return solveRead(currentGenes, activeIndex, seedLen, reads[i], fillEditDistance);
          }));
    } 

    for (int i = 0; i < reads.size(); ++i) {
      results[i].wait();
    }
    fprintf(stderr, "done (%.2lfs)\n", (clock() - starting_time) / double(CLOCKS_PER_SEC));
  }
}

void solveReads(Database& db, 
                const int seedLen,
                vector<shared_ptr<Read> >& reads) {
  solveReads(db, seedLen, reads, false);
  if (FLAGS_calc_edit_distance) {
    solveReads(db, seedLen, reads, true);
  }
}


void printStats(const vector<shared_ptr<Read> >& reads, const string& what) {
  map<int, int> stats;
  for (int i = 0; i < reads.size(); ++i) {
    shared_ptr<Read> read = reads[i];

    int mappingQuality = -1;
    if (what == "wgsim") {
      mappingQuality = read->validateWgsimMapping();
    } else {
      assert(0);
    }

    ++stats[mappingQuality];

    if (mappingQuality == 0) { // tocan je na prvom mjestu!
      double confidence1 = read->topMapping().score/read->size();
      double confidence2 = 1e10;
      
      // confidence2 trenutno nema toliku ulogu jer
      // na jednom genu/genomu mi podrazumijevamo samo jednu
      // poziciju
      //
      // TODO: razmatranje vise pozicija na genu/genomu moze
      // se rijesiti dodatnom podjelom svakog genoma/gena na
      // blokove pa u svakom bloku mozemo naci jedinstvenu poziciju
      // i uzeti blok s najvecim scoreom
      if (read->topMappings().size() >= 2) {
        OneMapping top1 = *(read->topMappings().begin());
        OneMapping top2 = *(++read->topMappings().begin());
        confidence2 = top1.score / top2.score;
      }

      if (confidence1 > FLAGS_confidence) {
        ++stats[-3];
      }
      if (confidence2 > 0.2) {
        ++stats[-4];
      }
    }

    bool validateFull = FLAGS_validate_wgsim == "full";
      
    if (mappingQuality != 0 && validateFull) {      
      printf("READ #%04d:\n", i);
      printf("mappingQuality = %d\n", mappingQuality);
      printf("id: %s\n", read->id().c_str());
      printf("data: %s\n", read->data().c_str());
      printf("mappings (%d):\n", (int)read->topMappings().size());
      
      for (auto mapping : read->topMappings()) {
        mapping.print(stdout);
      }
      puts("END READ\n");
    }
  }

  printf("Total reads: %d\n", (int)reads.size());
  printf("Number of reads not mapped: %d\n", stats[-1]);
  printf("Not even one kmer found: %d\n", stats[-2]);
  printf("Reads mapped with high confidence: %d\n", stats[-3]);
  printf("Reads mapped with high confidence2: %d\n", stats[-4]);

  for (map<int, int>::iterator it = stats.begin(); it != stats.end(); ++it) {
    if (it->first >= 0) {
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

    for (auto onemap : read->topMappings()) {
      fprintf(resultOut, ";%s,%lf,%d,%d,%d,%d", 
              onemap.geneDescriptor.c_str(),
      	      onemap.score, onemap.geneBegin, 
      	      onemap.geneEnd, onemap.isRC, 
              onemap.editDistance);
    }
    fprintf(resultOut, "\n");
  }

  fclose(resultOut);

  if (FLAGS_validate_wgsim != "") {
    printStats(reads, "wgsim");
  }
}

void clusterSolve(int argc, char* argv[]) {
  if (argc != 4) {
    printUsageAndExit();
  }

  string databasePath = argv[0];
  string indexPath = argv[1];
  string readsInputPath = argv[2];
  string readsOutputPath = argv[3];

  assert(isValidInputFile(databasePath));
  assert(isValidFolder(indexPath));
  assert(isValidInputFile(readsInputPath));
  assert(isValidOutputFile(readsOutputPath));

  MPI_Init(&argc, &argv);

  int noWorkers, myId;
  MPI_Comm_size(MPI_COMM_WORLD, &noWorkers);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);

  unsigned long long myChunkBegin = 1;
  unsigned long long myChunkEnd = 0;
  
  if (myId == 0) { // distributer
    // ovaj ce ucitavati s NFS-a i slati
    vector<unsigned long long> splitPos;
    splitReadInputFile(&splitPos, readsInputPath, noWorkers);

    assert((int)splitPos.size()-1 <= noWorkers);
    myChunkBegin = splitPos[0]; // i master ce odraditi dio
    myChunkEnd   = splitPos[1];

    for (int i = 1; i+1 < (int)splitPos.size(); ++i) {
      MPI_Send(&splitPos[i], 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Send(&splitPos[i+1], 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
    }
    for (int i = (int)splitPos.size()-1; i < noWorkers; ++i) {
      unsigned long long dummyBegin = 1, dummyEnd = 0;
      MPI_Send(&dummyBegin, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
      MPI_Send(&dummyEnd, 1, MPI_UNSIGNED_LONG_LONG, i, 0, MPI_COMM_WORLD);
    }
  } else {
    // ostali ce primiti
    MPI_Status stat;
    MPI_Recv(&myChunkBegin, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &stat);
    MPI_Recv(&myChunkEnd, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &stat);
  }
  
  // svatko ucita svoj dio ...
  vector<shared_ptr<Read> > reads;

  // ... odradi ...
  if (myChunkBegin < myChunkEnd) {
    inputReadsFileChunk(&reads, readsInputPath, myChunkBegin, myChunkEnd);
    
    Database db(databasePath, indexPath, FLAGS_seed_len, true);
    solveReads(db, FLAGS_seed_len, reads);
  }

  // ... i ispise u privremeni file
  int lastSlash = readsOutputPath.find_last_of("/");
  string tmpOutputFolder = "";
  if (lastSlash != string::npos) {
    tmpOutputFolder = readsOutputPath.substr(0, lastSlash) + "/";
  } 
  ostringstream tmpOutputFile;
  tmpOutputFile << tmpOutputFolder << myId;
  printf("tmpout: %s\n", tmpOutputFile.str().c_str()); fflush(stdout);
  printReads(reads, tmpOutputFile.str());
  
  MPI_Barrier(MPI_COMM_WORLD); // svi moraju dovrsiti ...

  // ... prije nego sto distributor pokupi sve privremene
  // fileove u konacni output cijelog clustera i obrise
  // privremene podatke
  if (myId == 0) {
    system(("rm -f " + readsOutputPath).c_str());
    for (int i = 0; i < noWorkers; ++i) {
      ostringstream ithFile;
      ithFile << tmpOutputFolder << i;
      system(("cat " + ithFile.str() + " >> " + readsOutputPath).c_str());
      system(("rm " + ithFile.str()).c_str());
    }
  }
  
  MPI_Finalize();
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
    } else {
      assert(isValidInputFile(argv[2]));
    }

    createIndex(argv[2], argv[3]);
  } else if (command == "test") {
    string readsFile = argv[2];
    
    vector<unsigned long long> filePos;
    splitReadInputFile(&filePos,
                       readsFile,
                       12);
    unsigned long long chunkChecksum = 0;
    int chunkTotalSize = 0;

    for (int i = 1; i < filePos.size(); ++i) {
      printf("[%llu %llu>\n", filePos[i-1], filePos[i]);

      vector<shared_ptr<Read> > chunkRead;
      inputReadsFileChunk(&chunkRead, readsFile,
                          filePos[i-1], filePos[i]);
      chunkTotalSize += chunkRead.size();
      for (int j = 0; j < chunkRead.size(); ++j) {
        chunkChecksum = chunkChecksum * 10007 + chunkRead[j]->checksum();
      }
    }
    printf("%d %llu\n", chunkTotalSize, chunkChecksum);

    vector<shared_ptr<Read> > allReads;
    inputReads(&allReads, readsFile);

    unsigned long long allChecksum = 0;
    for (int i = 0; i < allReads.size(); ++i) {
      allChecksum = allChecksum * 10007 + allReads[i]->checksum();
    }

    printf("%d %llu\n", (int)allReads.size(), allChecksum);
  } else if (command == "solve") {
    if (argc != 6) {
      printUsageAndExit();
    } else {
      assert(isValidInputFile(argv[2]));
      assert(isValidInputFile(argv[4]));
      assert(isValidOutputFile(argv[5]));
    }

    Database db(argv[2], argv[3], FLAGS_seed_len, true);
    vector<shared_ptr<Read> > reads;
    inputReads(&reads, argv[4], FLAGS_no_reads); 
    solveReads(db, FLAGS_seed_len, reads);
    printReads(reads, argv[5]);
  } else if (command == "cluster") {
    clusterSolve(argc-2, argv+2);
  } else {
    printUsageAndExit();
  }


  return 0;
}
