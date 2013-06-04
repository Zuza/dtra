// Class represents a single read to be aligned to a gene.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_READ
#define MAPPER_READ

#include <sys/types.h>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "core/bioinf_util.h"

struct OneMapping {
  double score;
  int genePos, isRC;

  int geneIdx;
  std::string geneDescriptor;
  std::string geneSegment; // optional

  OneMapping(double score, int genePos, int isRC, int geneIdx,
             std::string geneDescriptor, std::string geneSegment) :
    score(score), genePos(genePos), isRC(isRC), geneIdx(geneIdx),
    geneDescriptor(geneDescriptor), geneSegment(geneSegment) {}

  void print(FILE* out) {
    fprintf(out, "on gene %s (idx=%d), at position %d (RC=%d), score=%lf\n",
            geneDescriptor.c_str(), geneIdx, genePos, isRC, score);

    if (geneSegment.size() > 0) {
      fprintf(out, "segment: %s\n", geneSegment.c_str());
    }
  }

  bool operator < (const OneMapping& m) const {
    return score < m.score;
  }
};

class Read {
public:
  bool read(FILE* fi);
  void print(FILE* out);

  const std::string& id() const { return id_; }
  const std::string& data() const { return data_; }
  const int size() const { return data_.size(); }
  const OneMapping& topMapping() { return *topMappings_.rbegin(); }
  const std::multiset<OneMapping>& topMappings() { return topMappings_; }

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

  void updateMapping(double score, int genePos, int isRC, int geneIdx,
                     std::string geneDescriptor, 
                     std::string geneSegment = "");

  // ovo podrazumijeva da je ucitani read zapravo napravljen
  // simulatorom pa u liniji s imenom sadrzi i stvarne pozicije
  // na kojima se nalazi u genu
  int validateWgsimMapping();

  unsigned long long checksum() {
    unsigned long long ret = 0;
    for (int i = 0; i < id_.size(); ++i) ret = ret * 3137 + id_[i];
    for (int i = 0; i < data_.size(); ++i) ret = ret * 17 + data_[i];
    return ret;
  }

 private:
  std::string id_;
  std::string data_;
  
  //std::vector<OneMapping> topMappings_;
  std::multiset<OneMapping> topMappings_;
};

void splitReadInputFile(std::vector<unsigned long long>* filePos,
			const std::string& filePath,
			const int noParts);

void inputReadsFileChunk(std::vector<std::shared_ptr<Read> >* reads,
			 const std::string& filePath,
			 const unsigned long long& begin,
			 const unsigned long long& end);

#endif  // MAPPER_READ
