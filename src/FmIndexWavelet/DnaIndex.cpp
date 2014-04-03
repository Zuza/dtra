#include "DnaIndex.hpp"
#include "util.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

typedef unsigned long long ullint;

DnaIndex::DnaIndex(ullint all_genes_size) : all_genes_size(all_genes_size){
  index_created = false;
  data_store = new char[all_genes_size];
  data_ptr = 0;
  num_genes = 0;
}

DnaIndex::DnaIndex(FILE* in) {
  deserialize(in);
}

DnaIndex::~DnaIndex() {
  delete[] data_store; data_store = NULL;  
}

bool DnaIndex::insert_gene(const char* gene, ullint gene_len) {
  assert(!index_created);

  // using '!' as a gene separator
  // check that gene is over the alphabet
  for (ullint i = 0; i < gene_len; ++i) {
    assert(gene[i] == 'A' || gene[i] == 'C' || gene[i] == 'G' || gene[i] == 'T');
  }

  if (data_ptr + gene_len + 1 <= all_genes_size) {
    // this is ok, pass
  } else {
    // not enough memory
    return false;
  }

  ++num_genes;
  gene_starting_pos.push_back(data_ptr);
  // insert gene in data + separator '!'
  for (ullint i = 0; i < gene_len; ++i) {
    data_store[data_ptr++] = gene[i];
  }
  data_store[data_ptr++] = '!';
  assert(data_ptr <= all_genes_size);

  return true;
}

void DnaIndex::create_index() {
  const char alphabet[] = "ACGT!"; const int alphabet_sz = strlen(alphabet);

  // create the gene positions - before data_store is garbled
  alphabet_to_idxs(data_store, data_ptr, alphabet, alphabet_sz);
  
  separator_pos.reset(
    new RankedBitmap(
      data_store, data_ptr,
      /* everything is in the following value interval */ data_ptr,
      0, 4, 5 // "ACGT" = 0, '!' = 4
    ));

  idxs_to_alphabet(data_store, data_ptr, alphabet, alphabet_sz);

  // create index - careful! garbes data_store!
  fmindex.reset(new FmIndex(data_store, data_ptr, alphabet, alphabet_sz));

  index_created = true;
  delete[] data_store; data_store = NULL;
}

void DnaIndex::get_substring_pos(vector<pair<ullint, ullint> >& results, const char* query, int query_len, int limit) {
  assert(index_created);
  thread_local vector<ullint> raw_positions; raw_positions.clear(); // TODO: change to thread local!!!
  fmindex->get_substring_pos(raw_positions, query, query_len, limit);
  results.clear();
  for (ullint i = 0; i < raw_positions.size(); ++i) {
    ullint raw_pos = raw_positions[i];
    ullint gene_idx = raw_position_to_gene_idx(raw_pos);
    results.push_back(make_pair(gene_idx, raw_pos - gene_starting_pos[gene_idx]));
  }
}

ullint DnaIndex::raw_position_to_gene_idx(ullint raw_position) {
  return separator_pos->get_rank(1, raw_position);
}

void DnaIndex::serialize(FILE* out) const {
  assert(index_created);
  ::serialize(out, all_genes_size);
  ::serialize(out, num_genes);
  fmindex->serialize(out);
  serialize_vector(out, gene_starting_pos);
  separator_pos->serialize(out);
}

void DnaIndex::deserialize(FILE* in) {
  ::deserialize(in, all_genes_size);
  data_store = NULL;
  data_ptr = 0;
  index_created = true;
  ::deserialize(in, num_genes);
  fmindex.reset(new FmIndex(in));
  deserialize_vector(in, gene_starting_pos);
  separator_pos.reset(new RankedBitmap(in));
}
