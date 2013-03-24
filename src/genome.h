// Simple class that represents a genome. Contains the data, name and
// functions to read it from FASTA file.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_GENOME
#define MAPPER_GENOME

#include <cstring>
#include <string>
#include <vector>

class Genome {
 public:

  const std::string& name() const { return name_; }
  const std::string& data() const { return data_; }
  const int size() const { return data_.size(); }
  
  //private:
  std::string name_;
  std::string data_;
};

bool readGenome(Genome* g, FILE* inputFilePointer);

#endif  // MAPPER_GENOME
