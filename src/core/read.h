// Class represents a single read to be aligned to a gene.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_READ
#define MAPPER_READ

#include <sys/types.h>
#include <memory>
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

  unsigned long long checksum() {
    unsigned long long ret = 0;
    for (int i = 0; i < id_.size(); ++i) ret = ret * 3137 + id_[i];
    for (int i = 0; i < data_.size(); ++i) ret = ret * 17 + data_[i];
    return ret;
  }

 private:
  std::string id_;
  std::string data_;
  
  std::vector<OneMapping> topMappings_;
};

void splitReadInputFile(std::vector<unsigned long long>* filePos,
			const std::string& filePath,
			const int noParts);

void inputReadsFileChunk(std::vector<std::shared_ptr<Read> >* reads,
			 const std::string& filePath,
			 const unsigned long long& begin,
			 const unsigned long long& end);

#endif  // MAPPER_READ
