#include "core/database.h"
#include <cstdlib>
using namespace std;

const size_t kMaxBlockSize = 10000000; // 10 MB

Database::Database(const string& indexFolderPath,
                   const int seedLen) : seedLen_(seedLen),
                                        indexFolderPath_(indexFolderPath)
{
  assert(indexFolderPath.back() == '/');  
  dbFilePointer_ = NULL;
}


Database::Database(const string& databasePath,
                   const string& indexFolderPath,
                   const int seedLen) : seedLen_(seedLen),
                                        indexFolderPath_(indexFolderPath)                                        
{
  assert(indexFolderPath.back() == '/');

  dbFilePointer_ = fopen(databasePath.c_str(), "rt");
  assert(dbFilePointer_);

  currentIndex_ = 0;
}

Database::~Database() {
  if (dbFilePointer_) {
    fclose(dbFilePointer_);
  }
}

bool Database::readDbStoreIndex() { 
  clear_statistics();
  currentBlockNoBytes_ = 0;

  shared_ptr<Index> in(new Index(seedLen_));

  for (int iter = 0; ; ++iter) {
    shared_ptr<Gene> g(new Gene());
    if (!readGene(g.get(), dbFilePointer_)) {
      break;
    }

    in->insertGene(g.get());
    update_statistics(g.get());

    currentBlockNoBytes_ += g->name_.size();
    currentBlockNoBytes_ += g->data_.size();
    if (currentBlockNoBytes_ > kMaxBlockSize) {
      break;
    }
  }

  if (currentBlockNoBytes_ == 0) {
    char filename[100]; sprintf(filename, "%scount.txt", indexFolderPath_.c_str()); 
    FILE* out = fopen(filename, "w"); // write out how many index files there are
    assert(out);
    fprintf(out, "%d", currentIndex_);
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
  
  return true;
}

shared_ptr<Index> Database::readIndexFile(int which) {
  char filename[100]; sprintf(filename, "%s%d", indexFolderPath_.c_str(), which); 
  FILE* indexFile = fopen(filename, "rb");
  BufferedBinaryReader reader(indexFile);

  shared_ptr<Index> ptr(new Index(seedLen_));
  ptr->readIndex(reader);
  
  fclose(indexFile);
  assert(indexFile);
  return ptr;
}

int Database::getIndexFilesCount() {
  char filename[100]; sprintf(filename, "%scount.txt", indexFolderPath_.c_str()); 
  FILE* in = fopen(filename, "r"); // write out how many index files there are
  assert(in);
  int count; fscanf(in, "%d", &count);
  fclose(in);
  return count;
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
