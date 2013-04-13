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

class Index {
public:
  // Create index for a set of genes. 'seed_length' determines the
  // indexed unit length (substring).

  Index(int seedLength);

  void insertGene(Gene* gene);
  void prepareIndex(); // sort the hashes

  void getPositions(std::vector<std::pair<unsigned int, unsigned int> >* retVal, unsigned int hash); // hash_t

  // serialize & deserialize
  void writeIndex(BufferedBinaryWriter& writer);
  void readIndex(BufferedBinaryReader& reader);

  //  unsigned long long checksum();

private:
  std::pair<unsigned int, unsigned int> position_to_gene_position(unsigned int position);

  struct Entry {
    unsigned int position;
    unsigned int hash; // hash_t

    friend bool operator < (const Entry& a, const Entry& b) {
      return a.hash < b.hash;
    }
  };

  std::vector<int> geneStartingPos_;
  std::vector<Index::Entry> index_;
  int startingPos_;
  int seedLength_;
  bool indexPrepared_;
};


#endif
