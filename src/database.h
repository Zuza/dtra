#ifndef INNOCENTIVE_DATABASE
#define INNOCENTIVE_DATABASE

#include <string>
#include <vector>
#include <memory>

#include "genome.h"
#include "index.h"

class Database {
 public:
  /**
   * @param createIndex denotes whether index should be
   *                    read from previously made index file
   *                    or the purpose of this instance is
   *                    to create index and dump it into a file
   *                    (in both cases location is in indexFilePath)
   */
  Database(const std::string& databasePath,
	   const std::string& indexFilePath,
	   const bool createIndex);
  ~Database();

  bool readNextBlock();
  void printNames();
  size_t getCurrentBlockNoBytes() {
    return currentBlockNoBytes_;
  }

  unsigned long long checksum() {
    unsigned long long ret = 0;
    for (int i = 0; i < speciesIndex_.size(); ++i) {
      ret = ret * 1000000007 + speciesIndex_[i]->checksum();
    }
    return ret;
  }

 private:

 private:
  std::vector<std::shared_ptr<Genome> > species_;
  std::vector<std::shared_ptr<Index> > speciesIndex_;

  FILE* dbFilePointer_;
  FILE* indexFilePointer_;
  size_t currentBlockNoBytes_;
  bool createIndex_;
};

#endif
