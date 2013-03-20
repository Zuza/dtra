#ifndef INNOCENTIVE_DATABASE
#define INNOCENTIVE_DATABASE

#include <string>
#include <vector>

#include "genome.h"

class Database {
 public:
  Database(const std::string& databasePath);
  ~Database();

  bool readNextBlock();
  
 private:

 private:
  std::vector<Genome> species_;
  FILE* inputFilePointer_;
};

#endif
