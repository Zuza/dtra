#ifndef _DNA_INDEX_H_
#define _DNA_INDEX_H_

#include "FmIndex.hpp"
#include "RankedBitmap.hpp"

#include <memory>
#include <unordered_map>

typedef unsigned long long ullint;

class DnaIndex {
public:
  /*!
   * @param all_genes_size Cumulative size of all genes (without any metadata) in bytes
   */
  DnaIndex(ullint all_genes_size);

  // read index from FILE*
  DnaIndex(FILE* in);
  ~DnaIndex();

  // insert a new gene in the library in order
  // returns 1 on success, 0 on not enough memory
  // gene should be over the alphabet {A, C, G, T}
  bool insert_gene(const char* gene, ullint gene_len);

  // call after inserting all the genes (no more are allowed)
  // must be called before any substring query
  void create_index();

  // ----- for functions below this line, create_index() must be called beforehand ----- //

  // returns the vector of tuples (gene number, position in gene) of each location
  // where the query occurs in the gene library
  // at most limit results are returned
  // the results do NOT have to be sorted in any way
  // clears results
  void get_substring_pos(std::vector<pair<ullint, ullint> >& results, const char* query, int query_len, int limit = 100);

  void serialize(FILE* out) const;
  void deserialize(FILE* in);

private:

  ullint raw_position_to_gene_idx(ullint raw_position);

  ullint all_genes_size;
  char* data_store;  // used only before construction
  ullint data_ptr;   // used only before construction
  bool index_created;
  ullint num_genes;
  std::shared_ptr<FmIndex> fmindex;
  std::vector<ullint> gene_starting_pos;
  std::shared_ptr<RankedBitmap> separator_pos;  // separator positions, so we can map position -> (gene num, position in gene)
};

#endif // !_DNA_INDEX_H_
