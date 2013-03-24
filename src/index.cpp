#include "index.h"

using namespace std;

Index::Index(Genome* genome, int seed_length, int max_seed_positions)
  : seed_length_(seed_length), max_seed_positions_(max_seed_positions) {
  CreateIndex(genome);
}

void Index::CreateIndex(Genome* genome) {
  index_.clear();
  for (int i = 0; i <= genome->size() - seed_length_; ++i) {
    index_.insert(std::make_pair(genome->data().substr(i, seed_length_), i));
  }
}

bool Index::Query(const string& seed, vector<int>* res) {
  auto eq = index_.equal_range(seed);
  int count = distance(eq.first, eq.second);
  
  if (count > max_seed_positions_ || count <= 0) {
    return false;
  }
  
  for (auto it = eq.first; it != eq.second; ++it) {
    res->push_back(it->second);
  }
  
  return true;
}
