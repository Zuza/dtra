#include "core/database.h"
#include "FmIndexWavelet/DnaIndex.hpp"

#include <gflags/gflags.h>
#include <cstdlib>
#include <ctime>

using namespace std;

DEFINE_int32(indexPartSize, 100, "Sequence size in MBs that will go to"
             "one index file (only used during index construction)");

Database::Database(const string& databasePath, // ex. FASTA file
                   const string& indexFolderPath,
                   const int seedLen,
                   const bool indexCreated) : seedLen_(seedLen),
                                              indexFolderPath_(indexFolderPath),
                                              indexCreated_(indexCreated) {
  if (indexFolderPath_.back() != '/') {
    indexFolderPath_ += "/";
  }

  system(("mkdir -p " + indexFolderPath_).c_str());

  dbFilePointer_ = fopen(databasePath.c_str(), "rt");
  assert(dbFilePointer_);

  currentIndex_ = 0;
  if (indexCreated) {
    read_index_summaries();
  }
}

Database::~Database() {
  fclose(dbFilePointer_);
}

void Database::write_one_block(const DnaIndex& dna_index) {
  assert(indexFolderPath_.size() < 90);
  char filename[100]; sprintf(filename, "%s%d", indexFolderPath_.c_str(), currentIndex_);

  FILE* indexFile = fopen(filename, "wb");
  assert(indexFile);
  dna_index.serialize(indexFile);
  fclose(indexFile);
}

void Database::readDbStoreIndex() { 
  assert(indexCreated_ == false);
  currentBlockNoBytes_ = 0;
  clear_statistics();

  clock_t starting_time = clock();
  int num_genes = 0;
  long int starting_ftell = ftell(dbFilePointer_);
  shared_ptr<DnaIndex> dna_index;
  int last_percentage = 1;
  const size_t kMaxBlockSize = 1000000ull * FLAGS_indexPartSize;

  auto finish_one_block = [&]() {
    if (!dna_index) return;

    printf("time to read block = %.2lfs\n", double(clock() - starting_time) / CLOCKS_PER_SEC);
    starting_time = clock();

    dna_index->create_index();
    write_one_block(*dna_index);
    indexSummaries_.push_back(make_pair(starting_ftell, num_genes));

    printf("time to process block = %.2lfs\n", double(clock() - starting_time) / CLOCKS_PER_SEC);
    printf("longest entry = %llu bp, smallest entry = %llu bp, avg entry = %.2lf bp\n", maxSize_, minSize_, sizeSum_ / double(numGenes_));

    clear_statistics();
    starting_time = clock();
    num_genes = 0;
    starting_ftell = ftell(dbFilePointer_);
    last_percentage = 1;
    currentBlockNoBytes_ = 0;
    ++currentIndex_;
    dna_index.reset();    
  };

  do {
    shared_ptr<Gene> g(new Gene());
    if (!readGene(g.get(), dbFilePointer_)) {
      break;
    }
    g->subNonAcgtWithRandom(); // TODO: make this better?
    ++num_genes;

    do {
      if (!dna_index) {
        dna_index.reset(new DnaIndex(kMaxBlockSize));
      }
      bool success = dna_index->insert_gene(g->data(), g->dataSize());
      if (success) {
        update_statistics(g.get());
        currentBlockNoBytes_ += g->dataSize() + 1;
        break;
      } else {
        // this block is too small, add another
        finish_one_block();
      }
    } while(true);

    // print out completion percentage
    while (currentBlockNoBytes_*100 > last_percentage*kMaxBlockSize) {
      printf("%d%%.. ", last_percentage); fflush(stdout);
      ++last_percentage;
    }
  } while(true);
  finish_one_block();

  // write out the summary
  char filename[100]; sprintf(filename, "%s%s", indexFolderPath_.c_str(), "count.txt"); 
  FILE* out = fopen(filename, "wt"); // write out how many index files there are
  assert(out);

  fprintf(out, "%d\n", currentIndex_);
  for (auto p : indexSummaries_) {
    fprintf(out, "%ld %d\n", p.first, p.second);
  }
  fclose(out);
}

shared_ptr<DnaIndex> Database::readIndexFile(int which) {
  assert(indexCreated_ == true);
  char filename[100]; sprintf(filename, "%s%d", indexFolderPath_.c_str(), which); 

  FILE* indexFile = fopen(filename, "rb");
  assert(indexFile);
  shared_ptr<DnaIndex> ptr(new DnaIndex(indexFile)); // read it from file
  fclose(indexFile);

  clear_statistics();
  fseek(dbFilePointer_, indexSummaries_[which].first, SEEK_SET); // .first -> ftell of the index part
  genes_.clear();
  genes_.reserve(indexSummaries_[which].second); // .second -> number of genes in index part
  for (int i = 0; i < indexSummaries_[which].second; ++i) {
    shared_ptr<Gene> g(new Gene());
    assert(readGene(g.get(), dbFilePointer_));
    genes_.push_back(g);
    update_statistics(g.get());
  }

  return ptr;
}

int Database::getIndexFilesCount() {
  assert(indexCreated_ == true);
  return indexFilesCount_;
}

void Database::clear_statistics() {
  sizeSum_ = 0;
  minSize_ = 0-1; // OVERFLOW, MAX VALUE
  maxSize_ = 0;
  numGenes_ = 0;
}

void Database::update_statistics(Gene* gene) {
  size_t sz = gene->dataSize();
  minSize_ = min<unsigned long long>(minSize_, sz);
  maxSize_ = max<unsigned long long>(maxSize_, sz);
  sizeSum_ += sz;
  ++numGenes_;
}

void Database::read_index_summaries() {
  char filename[100]; sprintf(filename, "%scount.txt", indexFolderPath_.c_str());
  FILE* in = fopen(filename, "rt"); // write out how many index files there are

  if (!in) {
    fprintf(stderr, "FAILED reading number of index files from %s\n", filename);
    exit(1);
  }

  int count; fscanf(in, "%d", &count);
  indexFilesCount_ = count;
  indexSummaries_.resize(count);
  for (int i = 0; i < count; ++i) {
    fscanf(in, "%ld %d",
           &indexSummaries_[i].first, // ftell of the position
           &indexSummaries_[i].second); // number of genes
  }
  fclose(in);
}
