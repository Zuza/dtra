// Class is used to index and query the genome according to how SNAP works.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_INDEX
#define MAPPER_INDEX

#include <utility>
#include <tr1/unordered_map>

#include "genome.h"

class Index {
public:
  // Creates index for a given genome. 'seed_length' determines the
  // indexed unit length (substring) and 'max_seed_positions' the
  // number of seed hits beyond which results are ignored.
  Index(Genome* genome, int seed_length, int max_seed_positions);

  // Finds all positions in the genome for the given seed. If the
  // number of results exceeds 'max_seed_positions' or is 0, function
  // returns 'false'. Otherwise it returns 'true' with all positions
  // in 'ret'.  
  bool Query(const std::string& seed, std::vector<int>* res);
private:
  void CreateIndex(Genome* genome);

  std::tr1::unordered_multimap<std::string, int> index_;

  int seed_length_;
  int max_seed_positions_;
};


#endif
