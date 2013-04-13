#ifndef INNOCENTIVE_DATABASE
#define INNOCENTIVE_DATABASE

#include <string>
#include <vector>
#include <memory>

#include "BufferedBinaryReader.h"
#include "genome.h"
#include "index.h"

class Database {
 public:

  // read from database and store in index
  Database(const std::string& databasePath, // pass anything if you're only reading index files
           const std::string& indexFilePath,
           const int seedLen);

  // read from previously created index
  Database(const std::string& indexFolderPath,
           const int seedLen);

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

  /* unsigned long long checksum() { */
  /*   unsigned long long ret = 0; */
  /*   for (int i = 0; i < speciesIndex_.size(); ++i) { */
  /*     ret = ret * 1000000007 + speciesIndex_[i]->checksum(); */
  /*   } */
  /*   return ret; */
  /* } */

  const int& getSeedLen() const {
    return seedLen_;
  }

 private:
  void update_statistics(Genome* genome);
  void clear_statistics();

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

  unsigned long long sizeSum_;
  unsigned long long minSize_, maxSize_;
  unsigned long long numGenes_;
};

#endif
