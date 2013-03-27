// Class is used to index and query the genome according to how SNAP works.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_INDEX
#define MAPPER_INDEX

#include <utility>
#include <tr1/unordered_map>

#include "genome.h"

typedef unsigned long long keyType;
typedef unsigned long long indexType; // mozda ce i int biti dovoljan

class Index {
public:
  // Creates index for a given genome. 'seed_length' determines the
  // indexed unit length (substring) and 'max_seed_positions' the
  // number of seed hits beyond which results are ignored.
  Index(int seedLength);

  void create(Genome* genome);
  
  void appendToBinaryFile(FILE* fileToAppend);
  
  void readNextFromBinaryFile(FILE* indexInputFile);

  unsigned long long checksum() {
    unsigned long long ret = seedLength_;
    for (int i = 0; i < (int)positions_.size(); ++i) {
      ret = ret * 10007 + positions_[i].first;
      ret = ret * 3137  + positions_[i].second;
    }
    return ret;
  }
private:
  std::vector<std::pair<keyType, indexType> > positions_;

  int seedLength_;
};


#endif
