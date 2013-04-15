// Simple class that represents a gene. Contains the data, name and
// functions to read it from FASTA file.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_GENOME
#define MAPPER_GENOME

#include <cstring>
#include <string>
#include <vector>

class Gene {
 public:

  const std::string& name() const { return name_; }
  const char& name(size_t i) const { return name_[i]; }
  const std::string& data() const { return data_; }
  const char& data(size_t i) const { return data_[i]; }
  const std::string data(size_t a, size_t b) { 
    return std::string(data_.begin()+a, data_.begin()+b);
  }							
  const size_t size() const { return data_.size(); }

  void clear() { 
    name_.clear(); 
    data_.clear(); 
  }

  //private:
  std::string name_;
  std::string data_;
};

bool readGene(Gene* g, FILE* inputFilePointer);
bool printGene(Gene* g, FILE* outputFilePointer = stdout, int width = 80);

#endif  // MAPPER_GENOME
