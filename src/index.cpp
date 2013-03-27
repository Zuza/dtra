#include "index.h"
#include "bioinf_util.h"
#include <ctime>
#include <cstdlib>

#include <algorithm>
using namespace std;

namespace {
  unsigned leftPart(unsigned long long hsh) {
    return hsh >> 32;
  }
  unsigned rightPart(unsigned long long hsh) {
    return hsh & ((1LL<<32)-1);
  }
}

Index::Index(int seedLength) : seedLength_(seedLength) {
  assert(seedLength_ <= 32); // because we want the hash to fit 64bit integer
}

unsigned long long Index::checksum() {
  unsigned long long ret = seedLength_;
  for (int i = 0; i < (int)leftHashRanges_.size(); ++i) {
    ret = ret * 3137 + leftHashRanges_[i].first;
    ret = ret * 3137 + leftHashRanges_[i].second.first;
    ret = ret * 3137 + leftHashRanges_[i].second.second;
  }
  for (int i = 0; i < (int)rightHash_.size(); ++i) {
    ret = ret * 17 + rightHash_[i];
  }
  return ret;
}

void Index::create(Genome* genome) { 
  vector<pair<unsigned long long, int> > positions;

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
      positions.push_back(make_pair(rollingHash,
				    i+1-seedLength_));
    }
  }

  sort(positions.begin(), positions.end());

  leftHashRanges_.clear();
  rightHash_.clear();
  positions_.clear();

  if (positions.empty()) {
    // toliko je kratak da ni nema seedova
    return;
  }

  int first = 0;
  for (int i = 0; ; ++i) {
    if (i == (int)positions.size() ||
	leftPart(positions[first].first) != leftPart(positions[i].first)) {
      //assert(first < positions.size());
      leftHashRanges_.push_back(
      	make_pair(
      	  leftPart(positions[first].first),
      	  make_pair(first, i-1)
      	)
      );
      first = i;
      if (first == (int)positions.size()) {
	break;
      }
    }
    //assert(i<positions.size());
    rightHash_.push_back(rightPart(positions[i].first));
    positions_.push_back(positions[i].second);
  }
}

// TODO: buffered output
void Index::appendToBinaryFile(FILE* out) { 
  int leftSize = leftHashRanges_.size();
  fwrite(&leftSize, 4, 1, out);
  for (int i = 0; i < leftSize; ++i) {
    fwrite(&leftHashRanges_[i].first, 4, 1, out);
    fwrite(&leftHashRanges_[i].second.first, 4, 1, out);
    fwrite(&leftHashRanges_[i].second.second, 4, 1, out);
  }
  
  int rightSize = rightHash_.size();
  fwrite(&rightSize, 4, 1, out);
  for (int i = 0; i < rightSize; ++i) {
    fwrite(&rightHash_[i], 4, 1, out);
  }
  for (int i = 0; i < rightSize; ++i) {
    fwrite(&positions_[i], 4, 1, out);
  }
}

// TODO: buffered input
void Index::readNextFromBinaryFile(FILE* in) {
  int leftSize;
  fread(&leftSize, 4, 1, in);
  leftHashRanges_.resize(leftSize);

  for (int i = 0; i < leftSize; ++i) {
    fread(&leftHashRanges_[i].first, 4, 1, in);
    fread(&leftHashRanges_[i].second.first, 4, 1, in);
    fread(&leftHashRanges_[i].second.second, 4, 1, in);
  }

  int rightSize; 
  fread(&rightSize, 4, 1, in);
  rightHash_.resize(rightSize);
  positions_.resize(rightSize);

  for (int i = 0; i < rightSize; ++i) {
    fread(&rightHash_[i], 4, 1, in);
  }
  for (int i = 0; i < rightSize; ++i) {
    fread(&positions_[i], 4, 1, in);
  }
}


