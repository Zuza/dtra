#include <cassert>
#include <cstring>

#include <fstream>
#include <string>
#include <vector>

#include "read.h"
#include "util.h"
using namespace std;

DEFINE_int32(no_top_mappings, 5, "Number of best mappings a read records.");


bool Read::read(FILE* fi) {
  static char id[100001];
  static char data[100001];
  static char nesto[100001]; // TODO: sta je treca linija, za sad ignore
  static char kvaliteta[100001]; // TODO: cetvrta valjda kvalitetu oznacava, 
                                 // ignore za sada

  id[100000] = 0;
  if (fscanf(fi, "%s", id) != 1) return 0;
  assert(id[100000] == 0);

  data[100000] = 0;
  if(fscanf(fi, "%s", data) != 1) {
    printf("%s\n", id);
    assert(0);
  }
  assert(data[100000] == 0);

  nesto[100000] = 0;
  assert(fscanf(fi, "%s", nesto) == 1);
  assert(nesto[100000] == 0);

  kvaliteta[100000] = 0;
  assert(fscanf(fi, "%s", kvaliteta) == 1);
  assert(kvaliteta[100000] == 0);
  
  this->id_ = id;
  this->data_ = data;
  
  return 1;
}

void Read::print(FILE* out) {
  fprintf(out, "id: %s\n", id_.c_str());
  fprintf(out, "read: %s\n", data_.c_str());

  fprintf(out, "mappings:\n");
  for (int i = 0; i < topMappings_.size(); ++i) {
    topMappings_[i].print(out);
  }
}

int Read::validateFluxMapping(int maxOffset) {
  vector<string> tokens = Split(id_, ':');
  int pos1 = atoi(tokens[5].c_str()), pos2 = atoi(tokens[6].c_str());

  int geneIdx = atoi(tokens[0].c_str() + 1);
  --pos1; --pos2; // is this true?
  --geneIdx;
  
  for (int i = 0; i < topMappings_.size(); ++i) {
    bool geneMatch = false;

    if (geneIdx == topMappings_[i].geneIdx)
      geneMatch = true;

    if (geneMatch &&
        (abs(topMappings_[i].genePos-pos1) < maxOffset ||
         abs(topMappings_[i].genePos+size()-pos2) < maxOffset)) {
      return i;
    }
  }
  
  return -1;
}


int Read::validateWgsimMapping(int maxOffset) {
  vector<string> tokens = Split(id_, '|');
  int pos1 = -1000000, pos2 = -1000000;
  
  tokens[4] = tokens[4].substr(tokens[4].find("_"));
  assert(sscanf(tokens[4].c_str(), "_%d_%d", &pos1, &pos2) == 2);
  
  --pos1; --pos2;
  string readInGene = tokens[3];
  
  for (int i = 0; i < topMappings_.size(); ++i) {
    bool geneMatch = false;

    // ovo nije nuzno savrseno tocno, ali mislim da se u stvarnosti ne
    // dogadja slucaj kada ne radi
    if (topMappings_[i].geneDescriptor.find(readInGene) != string::npos) {
      geneMatch = true;
    }

    if (geneMatch &&
        (abs(topMappings_[i].genePos-pos1) < maxOffset ||
         abs(topMappings_[i].genePos-(pos2-size())) < maxOffset)) {
      return i;
    }
  }
  
  return -1;
}

void Read::updateMapping(double score, int genePos, int isRC, int geneIdx,
			 string geneDescriptor, 
			 string geneSegment) {
  topMappings_.push_back(OneMapping(score, genePos, isRC, geneIdx,
				    geneDescriptor, geneSegment));
  size_t i = topMappings_.size()-1;
  
  for ( ; i >= 1 && topMappings_[i-1] < topMappings_[i]; --i) {
    swap(topMappings_[i-1], topMappings_[i]);
  }
  
  if (topMappings_.size() > FLAGS_no_top_mappings) {
      topMappings_.pop_back();
  }
}


// TODO: mozda napraviti 'tezinsko' dijeljenje jer neke su 
// radilice jace od ostalih
void splitReadInputFile(vector<unsigned long long>* filePos,
			const string& filePath,
			const int noParts) {
  unsigned long long fileSize = getFileSize(filePath);
  unsigned long long chunkSize = fileSize / noParts;

  assert(filePos->empty());
  for (int i = 0; i < noParts; ++i) {
    filePos->push_back(i*chunkSize);
  }

  FILE* fi = fopen(filePath.c_str(), "rt");

  for (int i = 0; i < noParts; ++i) {
    assert(fseek(fi, (*filePos)[i], SEEK_SET) == 0);
    
    char previous, current;
    for (previous = '\n'; ; previous = current) {
      current = getc(fi);
      if (current == EOF) { 
	// ftell ce u ovom trenutku biti jednak 
	// velicini filea
	break;
      }

      if (current == '@' && previous == '\n') { 
	// pocetak nekog reada
	ungetc('@', fi);
	break;
      }
    }

    (*filePos)[i] = ftell(fi);
  }

  fclose(fi);

  filePos->push_back(fileSize);
  assert(is_sorted(filePos->begin(), filePos->end()));
  filePos->erase(unique(filePos->begin(), filePos->end()), filePos->end());

  assert((*filePos)[0] == 0);
  assert(filePos->back() == fileSize);
}

// ucitavam segment file izmedju [begin, end>
void inputReadsFileChunk(vector<shared_ptr<Read> >* reads,
			 const string& filePath,
			 const unsigned long long& begin,
			 const unsigned long long& end) {
  FILE* fi = fopen(filePath.c_str(), "rt");
  assert(fi);

  fseek(fi, begin, SEEK_SET);
  while (true) {
    //printf("%llu %llu\n", ftell(fi), end);
    shared_ptr<Read> tmp(new Read());

    if (!tmp->read(fi)) {
      break;
    }

    if (ftell(fi) <= end) {
      reads->push_back(tmp);
    } else {
      break;
    }
  }

  fclose(fi);
}
