/**
 * @authors: Goran Zuzic (zuza777@gmail.com)
 *           Filip Pavetic (fpavetic@gmail.com)
 */

#ifndef MAPPER_GENOME
#define MAPPER_GENOME

#include "core/bioinf_util.h"

#include <cstring>
#include <string>
#include <vector>

class Gene {
 public:
  Gene() { name_ = data_ = 0; }
  ~Gene() { clear(); }

  const char* name() const { return name_; }
  const char name(size_t i) const { return name_[i]; }
  const char* data() const { return data_; }
  const char data(size_t i) const { return data_[i]; }
  const std::string& description() const { return description_; }

  const size_t dataSize() const { return data_len_; }
  const size_t nameSize() const { return name_len_; }

  void subNonAcgtWithRandom() {
    ::subNonAcgtWithRandom(data_, data_len_);
  }

  void clear() { 
    if (name_) { 
      free(name_);
    }
    if (data_) {
      free(data_);
    }
    name_ = data_ = 0;
    name_len_ = 0;
    data_len_ = 0;
  }

  friend bool readGene(Gene* g, FILE* inputFilePointer);

 private:
  void createDescription();

 private:                                   
  char* name_; size_t name_len_;
  char* data_; size_t data_len_;
  std::string description_;
};

bool readGene(Gene* g, FILE* inputFilePointer);
size_t printGene(Gene* g, FILE* outputFilePointer = stdout, int width = 80);

#endif  // MAPPER_GENOME
