#ifndef INNOCENTIVE_DATABASE
#define INNOCENTIVE_DATABASE

#include <string>
#include <vector>
#include <memory>

#include "BufferedBinaryReader.h"
#include "gene.h"
#include "index.h"

DECLARE_bool(discard_freq_seeds);

class Database {
 public:

  Database(const std::string& databasePath, // pass anything if you're only reading index files
           const std::string& indexFilePath,
           const int seedLen,
           const bool indexCreated);

  ~Database();

  bool readDbStoreIndex();
  int getIndexFilesCount();
  std::shared_ptr<Index> readIndexFile(int which);

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

  int indexFilesCount_;

  FILE* dbFilePointer_;
  size_t currentBlockNoBytes_;
  bool createIndex_;
  int seedLen_;
  int currentIndex_;

  std::string indexFolderPath_;

  // koristi se kod izgradnje indeksa
  std::shared_ptr<BufferedBinaryWriter> bufferedWriter_;

  // koristi se kod ucitavanja indeksa
  std::shared_ptr<BufferedBinaryReader> bufferedReader_;
  
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
