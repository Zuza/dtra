// Class represents a single read to be aligned to a gene.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_READ
#define MAPPER_READ

#include <sys/types.h>
#include "core/bioinf_util.h"

#include <string>

const int kNoTopMappings = 10;

struct OneMapping {
  double score;
  int geneId, genePos, isRC;

  OneMapping(double score, int geneId, int genePos, int isRC) :
    score(score), geneId(geneId), genePos(genePos), isRC(isRC) {}

  void print() {
    printf("on gene %d, at position %d (RC=%d), score=%lf\n",
	   geneId, genePos, isRC, score);
  }

  bool operator < (const OneMapping& m) {
    return score < m.score;
  }
};

class Read {
public:
  bool read(FILE* fi);
  void print();

  const std::string& id() const { return id_; }
  const std::string& data() const { return data_; }
  const int size() const { return data_.size(); }
  const std::vector<OneMapping>& topMappings() { return topMappings_; }
  const OneMapping& topMapping(size_t i) { return topMappings_[i]; }

  // ako je reverseComplement == true, 
  // onda ce vratiti i-tu bazu reverse
  // complementa, inace vracam normalno i-tu
  // bazu reada
  char get(const size_t& i, const bool reverseComplement = false) { 
    if (!reverseComplement) {
      return data_[i]; 
    }
    return getBaseComplement(data_[data_.size()-1-i]);
  }

  void updateMapping(double score, int geneId, int genePos, int isRC) {
    topMappings_.push_back(OneMapping(score, geneId, genePos, isRC));
    size_t i = topMappings_.size()-1;
    
    for ( ; i >= 1 && topMappings_[i-1] < topMappings_[i]; --i) {
      std::swap(topMappings_[i-1], topMappings_[i]);
    }

    if (topMappings_.size() > kNoTopMappings) {
      topMappings_.pop_back();
    }
  }


private:
  std::string id_;
  std::string data_;

  std::vector<OneMapping> topMappings_;
};

#endif  // MAPPER_READ
