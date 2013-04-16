// Class represents a single read to be aligned to a gene.
//
// Authors: Matija Osrecki, Filip Pavetic

#ifndef MAPPER_READ
#define MAPPER_READ

#include <sys/types.h>
#include <string>

#include "core/bioinf_util.h"
#include "core/util.h"

const int kNoTopMappings = 5;

struct OneMapping {
  double score;
  int geneId, genePos, isRC;

  std::string geneStrId;
  std::string geneSegment; // optional

  OneMapping(double score, int geneId, int genePos, int isRC, 
	     std::string geneName, std::string geneSegment) :
  score(score), geneId(geneId), genePos(genePos), isRC(isRC),
    geneStrId(geneName), geneSegment(geneSegment) {
    std::vector<std::string> tokens = Split(geneName, '|');
    geneStrId = tokens[3];
  }

  void print() {
    printf("on gene %d, at position %d (RC=%d), score=%lf\n",
	   geneId, genePos, isRC, score);
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

  void updateMapping(double score, int geneId, int genePos, int isRC, 
		     std::string geneName, std::string geneSegment = "") {
    topMappings_.push_back(OneMapping(score, geneId, genePos, isRC, 
				      geneName, geneSegment));
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
  int validateMapping(int maxOffset = 5) {
    std::vector<std::string> tokens = Split(id_, '|');
    int pos1 = -1000000, pos2 = -1000000;
    sscanf(tokens[4].c_str(), "_%d_%d", &pos1, &pos2); 
    --pos1; --pos2;
    std::string geneId = tokens[3];

    for (int i = 0; i < topMappings_.size(); ++i) {
      int dokle = topMappings_[i].geneStrId.find(" ");
      std::string topMapId = topMappings_[i].geneStrId.substr(0, dokle);

      if (geneId == topMapId &&
	  (abs(topMappings_[i].genePos-pos1) < maxOffset ||
	   abs(topMappings_[i].genePos-(pos2-size())) < maxOffset)) {
	return i;
      }
    }

    return -1;
  }

 private:
  std::string id_;
  std::string data_;

  std::vector<OneMapping> topMappings_;
};

#endif  // MAPPER_READ
