#ifndef INNOCENTIVE_DATABASE
#define INNOCENTIVE_DATABASE

#include <string>
#include <vector>

#include "genome.h"
#include "index.h"

class Database {
 public:
  Database(const std::string& databasePath);
  ~Database();

  bool readNextBlock();
  void printNames();
  size_t getCurrentBlockNoBytes() {
    return currentBlockNoBytes_;
  }

 private:

 private:
  std::vector<Genome> species_;
  std::vector<Index> speciesIndex_;

  FILE* inputFilePointer_;
  size_t currentBlockNoBytes_;
};

#endif
