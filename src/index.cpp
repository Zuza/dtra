#include "index.h"
#include "bioinf_util.h"
#include <ctime>
#include <cstdlib>

#include <algorithm>
using namespace std;

Index::Index(int seedLength) : seedLength_(seedLength) {
  assert(seedLength_ <= 32); // because we want the hash to fit 64bit integer
}

void Index::create(Genome* genome) {
  positions_.clear();

  unsigned long long subtractPower = 1;
  for (int i = 0; i < seedLength_; ++i) {
    subtractPower *= 4;
  }

  unsigned long long rollingHash = 0;

  for (size_t i = 0; i < genome->size(); ++i) {
    int next = baseToInt(genome->data()[i]);

    rollingHash = rollingHash * 4 + next;
    if (i >= seedLength_) {
      rollingHash %= subtractPower;
    }

    if (i+1 >= seedLength_) {
      positions_.push_back(make_pair(rollingHash,
				     i+1-seedLength_));
    }
  }

  sort(positions_.begin(), positions_.end());
}

void Index::appendToBinaryFile(FILE* out) {
  int sz = positions_.size();
  fwrite(&sz, 4, 1, out);
  for (int i = 0; i < (int)positions_.size(); ++i) {
    fwrite(&positions_[i].first, 8, 1, out);
    fwrite(&positions_[i].second, 8, 1, out);
  }
}

void Index::readNextFromBinaryFile(FILE* in) {
  int sz; 
  fread(&sz, 4, 1, in);
  positions_.resize(sz);

  for (int i = 0; i < sz; ++i) {
    unsigned long long hsh;
    unsigned long long pos; 
    fread(&hsh, 8, 1, in);
    fread(&pos, 8, 1, in);
    positions_[i] = make_pair(hsh, pos);
  }
}


