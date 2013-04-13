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
  /**
   * @param createIndex denotes whether index should be
   *                    read from previously made index file
   *                    or the purpose of this instance is
   *                    to create index and dump it into a file
   *                    (in both cases location is in indexFilePath)
   */
  Database(const std::string& databasePath,
	   const std::string& indexFilePath,
	   const int seedLen,
	   const bool createIndex);
  ~Database();

  bool readNextBlock();
  void printNames();
  size_t getCurrentBlockNoBytes() {
    return currentBlockNoBytes_;
  }
  size_t getCurrentBlockNoGenes() {
    return species_.size();
  }
  double getAverageGeneLength() {
    double avg = 0;
    for (int i = 0; i < (int)species_.size(); ++i) {
      avg += species_[i]->data().size();
    }
    return avg / getCurrentBlockNoGenes();
  }
  void getMinMaxGeneLength(size_t* minLen, size_t* maxLen) {
    *minLen = 1000000000000LL;
    *maxLen = 0;
    for (int i = 0; i < (int)species_.size(); ++i) {
      size_t len = species_[i]->data().size();
      *minLen = std::min(*minLen, len);
      *maxLen = std::max(*maxLen, len);
    }
  }

  unsigned long long checksum() {
    unsigned long long ret = 0;
    for (int i = 0; i < speciesIndex_.size(); ++i) {
      ret = ret * 1000000007 + speciesIndex_[i]->checksum();
    }
    return ret;
  }

  const int& getSeedLen() const {
    return seedLen_;
  }

 private:

 private:
  std::vector<std::shared_ptr<Genome> > species_;
  std::vector<std::shared_ptr<Index> > speciesIndex_;

  FILE* dbFilePointer_;
  FILE* indexFilePointer_;
  size_t currentBlockNoBytes_;
  bool createIndex_;
  int seedLen_;

  // koristi se kod izgradnje indeksa
  std::shared_ptr<BufferedBinaryWriter> bufferedWriter_;

  // koristi se kod ucitavanja indeksa
  std::shared_ptr<BufferedBinaryReader> bufferedReader_;

};

#endif
