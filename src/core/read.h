// Class represents a single read to be aligned to a gene.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_READ
#define MAPPER_READ

#include <sys/types.h>
#include <string>

#include "core/bioinf_util.h"

const int kNoTopMappings = 5;
const int kMaxOffset = 5;

struct OneMapping {
  double score;
  int genePos, isRC;

  std::string geneDescriptor;
  std::string geneSegment; // optional

  OneMapping(double score, int genePos, int isRC, 
             std::string geneDescriptor, std::string geneSegment) :
    score(score), genePos(genePos), isRC(isRC),
    geneDescriptor(geneDescriptor), geneSegment(geneSegment) {}

  void print() {
    printf("on gene %s, at position %d (RC=%d), score=%lf\n",
           geneDescriptor.c_str(), genePos, isRC, score);
    printf("segment: %s\n", geneSegment.c_str());
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

  void updateMapping(double score, int genePos, int isRC, 
                     std::string geneDescriptor, 
                     std::string geneSegment = "") {
    topMappings_.push_back(OneMapping(score, genePos, isRC, 
                                      geneDescriptor, geneSegment));
    size_t i = topMappings_.size()-1;
    
    for ( ; i >= 1 && topMappings_[i-1] < topMappings_[i]; --i) {
      std::swap(topMappings_[i-1], topMappings_[i]);
    }

    if (topMappings_.size() > kNoTopMappings) {
      topMappings_.pop_back();
    }
  }

  // ovo podrazumijeva da je ucitani read zapravo napravljen
  // simulatorom pa u liniji s imenom sadrzi i stvarne pozicije
  // na kojima se nalazi u genu
  int validateWgsimMapping(int maxOffset = kMaxOffset);

private:
  std::string id_;
  std::string data_;

  std::vector<OneMapping> topMappings_;
};

#endif  // MAPPER_READ
