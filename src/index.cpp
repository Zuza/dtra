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

void Index::appendToBufferedBinaryWriter(BufferedBinaryWriter& writer) { 
  int leftSize = leftHashRanges_.size();

  writer.writeSigned32(leftSize);
  for (int i = 0; i < leftSize; ++i) {
    writer.writeUnsigned32(leftHashRanges_[i].first);
    writer.writeSigned32(leftHashRanges_[i].second.first);
    writer.writeSigned32(leftHashRanges_[i].second.second);
  }
  
  int rightSize = rightHash_.size();
  writer.writeSigned32(rightSize);
  for (int i = 0; i < rightSize; ++i) {
    writer.writeUnsigned32(rightHash_[i]);
  }
  for (int i = 0; i < rightSize; ++i) {
    writer.writeSigned32(positions_[i]);
  }
  
  writer.flushBuffer(true);
}

void Index::readNextFromBufferedReader(BufferedBinaryReader& reader) {
  int leftSize; reader.readSigned32(&leftSize);
  leftHashRanges_.resize(leftSize);
  
  for (int i = 0; i < leftSize; ++i) {
    reader.readUnsigned32(&leftHashRanges_[i].first);
    reader.readSigned32(&leftHashRanges_[i].second.first);
    reader.readSigned32(&leftHashRanges_[i].second.second);
  }

  int rightSize; reader.readSigned32(&rightSize);
  rightHash_.resize(rightSize);
  for (int i = 0; i < rightSize; ++i) {
    reader.readUnsigned32(&rightHash_[i]);
  }

  positions_.resize(rightSize);
  for (int i = 0; i < rightSize; ++i) {
    reader.readSigned32(&positions_[i]);
  }
}


