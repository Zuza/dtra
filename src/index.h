// Class is used to index and query the genome according to how SNAP works.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_INDEX
#define MAPPER_INDEX

#include <utility>
#include <tr1/unordered_map>

#include "genome.h"

typedef unsigned long long keyType;

// u nt bazi najdulji gen je dugacak 150M pa pozicije stanu u int
typedef int indexType;

class Index {
public:
  // Creates index for a given genome. 'seed_length' determines the
  // indexed unit length (substring) and 'max_seed_positions' the
  // number of seed hits beyond which results are ignored.
  Index(int seedLength);

  void create(Genome* genome);
  
  void appendToBinaryFile(FILE* fileToAppend);
  
  void readNextFromBinaryFile(FILE* indexInputFile);

  unsigned long long checksum();

private:
  //std::vector<std::pair<keyType, indexType> > positions_;
  std::vector<std::pair<unsigned, std::pair<int, int> > > leftHashRanges_;
  std::vector<unsigned> rightHash_;
  std::vector<int> positions_;

  int seedLength_;
};


#endif
