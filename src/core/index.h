// Class is used to index and query the gene according to how SNAP works.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_INDEX
#define MAPPER_INDEX

#include <vector>

#include "BufferedBinaryReader.h"
#include "BufferedBinaryWriter.h"
#include "gene.h"

// typedef unsigned long long keyType;

typedef unsigned long long hash_t;

class Index {
public:
  // Create index for a set of genes. 'seed_length' determines the
  // indexed unit length (substring).

  Index(int seedLength);

  void insertGene(Gene* gene);
  void prepareIndex(); // sort the hashes

  // in: hash value
  // out: retval
  //   -> vector of pairs -> first is the gene number
  //                         second is the position with the gene
  void getPositions(std::vector<std::pair<unsigned int, unsigned int> >* retVal, hash_t hash); // hash_t

  const int& getSeedLen() { return seedLength_; }

  // serialize & deserialize
  void writeIndex(BufferedBinaryWriter& writer);
  void readIndex(BufferedBinaryReader& reader);

private:
  std::pair<unsigned int, unsigned int> position_to_gene_position(size_t position);

  struct Entry {
    size_t position;
    hash_t hash; // hash_t

    friend bool operator < (const Entry& a, const Entry& b) {
      return a.hash < b.hash;
    }
  };

  std::vector<size_t> geneStartingPos_;
  std::vector<Index::Entry> index_;
  size_t startingPos_;
  int seedLength_;
  bool indexPrepared_;
};


#endif
