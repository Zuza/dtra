// Class represents a single read to be aligned to a genome.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_READ
#define MAPPER_READ

#include <sys/types.h>

#include <string>

class Read {
public:
  bool read(FILE* fi);
  void print();

  const std::string& id() const { return id_; }
  const std::string& data() const { return data_; }
  const int size() const { return data_.size(); }
  const char& operator [] (const size_t& i) { return data_[i]; }

private:
  std::string id_;
  std::string data_;
};

#endif  // MAPPER_READ
