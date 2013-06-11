
/**
 * @authors: Goran Zuzic (zuza777@gmail.com)
 */

#ifndef MAPPER_INDEX
#define MAPPER_INDEX

#include <vector>
#include <gflags/gflags.h>

#include "BufferedBinaryReader.h"
#include "BufferedBinaryWriter.h"
#include "gene.h"

typedef unsigned long long hash_t;

class Index {
public:
  // Create index for a set of genes. 'seed_length' determines the
  // indexed unit length (substring).

  Index(int seedLength);

  struct iterator {
    iterator(const Index* idx, 
	     const hash_t& hash, 
	     const int& querySeedLen);

    void reset();
    bool done();
    void advance();
    std::pair<unsigned int, unsigned int> get();

    /* void setStartingPos(size_t where); */

    const Index* idx_;
    hash_t hash_;
    size_t begin, end, curr, currStartingPos;
  };

  friend class iterator;

  void insertGene(Gene* gene);
  void prepareIndex(); // sort the hashes

  // in: hash value
  // out: retval
  //   -> vector of pairs -> first is the gene number
  //                         second is the position with the gene
  Index::iterator getPositions(const hash_t& hash, 
			       const int& querySeedLen);

  const int& getSeedLen() const { return seedLength_; }

  // serialize & deserialize
  void writeIndex(BufferedBinaryWriter& writer);
  void readIndex(BufferedBinaryReader& reader);

  // optimization
  void discardFrequentSeeds();

private:
  std::pair<unsigned int, unsigned int> position_to_gene_position(size_t position) const;

  struct Entry {
    size_t position;
    hash_t hash;

    friend bool operator < (const Entry& a, const Entry& b) {
      return a.hash < b.hash;

      // zbog toga sto se koristi stable_sort, oni s jednakim
      // hashem bit ce sortirani po poziciji
    }
  };

  struct VariableSeedLenCmp {
    VariableSeedLenCmp(int targetSeedLen, int indexSeedLen) {
      trimLen = indexSeedLen - targetSeedLen;
    }
    
    int trimLen;
    
    bool operator() (const Index::Entry& a, const Index::Entry& b) {
      return (a.hash >> (2*trimLen)) < (b.hash >> (2*trimLen));
    }
  };

  std::vector<size_t> geneStartingPos_;
  std::vector<Index::Entry> index_;
  size_t startingPos_;
  int seedLength_;
  bool indexPrepared_;
};


#endif
