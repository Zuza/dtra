#ifndef INNOCENTIVE_DATABASE
#define INNOCENTIVE_DATABASE

#include <string>
#include <vector>
#include <memory>

#include "FmIndexWavelet/DnaIndex.hpp"
#include "gene.h"

class Database {
 public:

  Database(const std::string& databasePath, // ignored if you're only reading index files
           const std::string& indexFilePath,
           const int seedLen,
           const bool indexCreated);

  ~Database();

  // read from the fasta file and store the corresponding index
  void readDbStoreIndex();
  int getIndexFilesCount();
  std::shared_ptr<DnaIndex> readIndexFile(int which);

  size_t getCurrentBlockNoBytes() {
    return currentBlockNoBytes_;
  }
  size_t getCurrentBlockNoGenes() {
    return numGenes_;
  }
  double getAverageGeneLength() {
    return double(sizeSum_) / numGenes_;
  }
  void getMinMaxGeneLength(size_t* minLen, size_t* maxLen) {
    *minLen = minSize_;
    *maxLen = maxSize_;
  }

  const int& getSeedLen() const {
    return seedLen_;
  }

  std::vector<std::shared_ptr<Gene> >& getGenes() {
    return genes_;
  }

 private:
  void update_statistics(Gene* gene);
  void clear_statistics();
  void read_index_summaries();
  void write_one_block(const DnaIndex& dna_index);

  int indexFilesCount_;

  FILE* dbFilePointer_;
  size_t currentBlockNoBytes_;
  bool createIndex_;
  int seedLen_;
  int currentIndex_;

  std::string indexFolderPath_;

  // first = ftell of the index start in the database
  // second = number of genes in the index file
  std::vector<std::pair<long int, int> > indexSummaries_;

  // gene vector
  std::vector<std::shared_ptr<Gene> > genes_;

  unsigned long long sizeSum_;
  unsigned long long minSize_, maxSize_;
  unsigned long long numGenes_;
  bool indexCreated_;
};

#endif
