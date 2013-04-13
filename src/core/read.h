// Class represents a single read to be aligned to a genome.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_READ
#define MAPPER_READ

#include <sys/types.h>

#include <string>

using namespace std;

class Read {
public:
  int readFromStdin();
  void print();

  const string& id() const { return id_; }
  const string& data() const { return data_; }
  const int size() const { return data_.size(); }
  const char& operator [] (const size_t& i) { return data_[i]; }

private:
  string id_;
  string data_;
};

#endif  // MAPPER_READ
