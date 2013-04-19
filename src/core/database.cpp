#include "core/database.h"
#include <cstdlib>
using namespace std;

const size_t kMaxBlockSize = 100000000; // 100 MB

Database::Database(const string& databasePath,
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

bool Database::readDbStoreIndex() { 
  assert(indexCreated_ == false);
  clear_statistics();
  currentBlockNoBytes_ = 0;

  int num_genes = 0;
  long int starting_ftell = ftell(dbFilePointer_);
  shared_ptr<Index> in(new Index(seedLen_));

  for (int iter = 0; ; ++iter) {
    shared_ptr<Gene> g(new Gene());
    if (!readGene(g.get(), dbFilePointer_)) {
      break;
    }

    ++num_genes;
    in->insertGene(g.get());
    update_statistics(g.get());

    currentBlockNoBytes_ += g->nameSize();
    currentBlockNoBytes_ += g->dataSize();
    if (currentBlockNoBytes_ > kMaxBlockSize) {
      break;
    }
  }

  if (currentBlockNoBytes_ == 0) { // write the summary
    char filename[100]; sprintf(filename, "%s%s", indexFolderPath_.c_str(), "count.txt"); 
    FILE* out = fopen(filename, "w"); // write out how many index files there are
    assert(out);

    fprintf(out, "%d\n", currentIndex_);
    for (auto p : indexSummaries_) {
      fprintf(out, "%ld %d\n", p.first, p.second);
    }
    
    fclose(out);
    return false;
  }

  in->prepareIndex();
  assert(indexFolderPath_.size() < 90);
  char filename[100]; sprintf(filename, "%s%d", indexFolderPath_.c_str(), currentIndex_); 
  ++currentIndex_;
  
  FILE* indexFile = fopen(filename, "wb");
  assert(indexFile);
  BufferedBinaryWriter writer(indexFile);
  in->writeIndex(writer);
  fclose(indexFile);

  indexSummaries_.push_back(make_pair(starting_ftell, num_genes));
  return true;
}

shared_ptr<Index> Database::readIndexFile(int which) {
  assert(indexCreated_ == true);
  char filename[100]; sprintf(filename, "%s%d", indexFolderPath_.c_str(), which); 
  FILE* indexFile = fopen(filename, "rb");
  assert(indexFile);

  BufferedBinaryReader reader(indexFile);
  shared_ptr<Index> ptr(new Index(seedLen_));
  ptr->readIndex(reader);

  fseek(dbFilePointer_, indexSummaries_[which].first, SEEK_SET);
  genes_.clear();
  genes_.reserve(indexSummaries_[which].second);
  for (int i = 0; i < indexSummaries_[which].second; ++i) {
    shared_ptr<Gene> g(new Gene());
    assert(readGene(g.get(), dbFilePointer_));
    genes_.push_back(g);
  }
  
  fclose(indexFile);
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
  size_t sz = gene->size();
  minSize_ = min<unsigned long long>(minSize_, sz);
  maxSize_ = max<unsigned long long>(maxSize_, sz);
  sizeSum_ += sz;
  ++numGenes_;
}

void Database::read_index_summaries() {
  char filename[100]; sprintf(filename, "%scount.txt", indexFolderPath_.c_str());
  FILE* in = fopen(filename, "r"); // write out how many index files there are

  if (!in) {
    fprintf(stderr, "FAILED reading number of index files from %s\n", filename);
    exit(1);
  }

  int count; fscanf(in, "%d", &count);
  indexFilesCount_ = count;
  indexSummaries_.resize(count);
  for (int i = 0; i < count; ++i) {
    fscanf(in, "%ld %d",
           &indexSummaries_[i].first,
           &indexSummaries_[i].second);
  }
  fclose(in);
}
