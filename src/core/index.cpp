#include "core/index.h"
#include "core/bioinf_util.h"
#include <ctime>
#include <cstdlib>
#include <algorithm>

#include "BufferedBinaryWriter.h"

using namespace std;

Index::Index(int seedLength) : seedLength_(seedLength) {
  indexPrepared_ = false;
  startingPos_ = 0;
  assert(seedLength_ <= 16); // currently unsigned int
  assert(seedLength_ <= 32); // maybe later: because we want the hash to fit 64bit integer ?
}

void Index::insertGene(Gene* gene) { 
  // TODO: prealociraj sve!

  unsigned long long subtractPower = 1;
  for (int i = 0; i < seedLength_; ++i) {
    subtractPower *= 4;
  }

  unsigned long long rollingHash = 0;

  for (size_t i = 0; i < gene->size(); ++i) {
    int next = baseToInt(gene->data()[i]);

    rollingHash = rollingHash * 4 + next;
    if (i >= seedLength_) {
      rollingHash %= subtractPower; // TODO(zuza): ANDanje umjesto moda
    }

    if (i+1 >= seedLength_) {
      Index::Entry tmp;
      tmp.position = startingPos_ + i+1 - seedLength_;
      tmp.hash = rollingHash;
      index_.push_back(tmp);
    }
  }

  geneStartingPos_.push_back(startingPos_);
  startingPos_ += gene->size();
}

pair<unsigned int, unsigned int> Index::position_to_gene_position(unsigned int position) {
  auto it = upper_bound(geneStartingPos_.begin(), geneStartingPos_.end(), position);
  assert(it != geneStartingPos_.begin());
  --it;
  return make_pair<unsigned int, unsigned int>(it - geneStartingPos_.begin(), position - *it);
}

void Index::getPositions(vector<pair<unsigned int, unsigned int> >* retVal, unsigned int hash) {
  assert(indexPrepared_);
  assert(retVal->empty());
  Entry tmp; tmp.hash = hash;
  auto pair_lb_ub = equal_range(index_.begin(), index_.end(), tmp);
  for ( ; pair_lb_ub.first != pair_lb_ub.second; ++pair_lb_ub.first) {
    Entry& entry = *pair_lb_ub.first;
    retVal->push_back(position_to_gene_position(entry.position));
  }
}

void Index::prepareIndex() {
  // positions are in increasing order
  stable_sort(index_.begin(), index_.end());
  indexPrepared_ = true;
}

void Index::writeIndex(BufferedBinaryWriter& writer) {
  assert(indexPrepared_);
  writer.writeSigned32(seedLength_);
  writer.writeVector(geneStartingPos_);
  writer.writeVector(index_);
}

void Index::readIndex(BufferedBinaryReader& reader) {
  indexPrepared_ = true;
  int storedSeedLength; reader.readSigned32(&storedSeedLength);
  assert(seedLength_ == storedSeedLength);
  reader.readVector(geneStartingPos_);
  reader.readVector(index_);
}
