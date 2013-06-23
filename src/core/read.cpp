#include <cassert>
#include <cstring>

#include <fstream>
#include <string>
#include <vector>

#include "read.h"
#include "util.h"
using namespace std;

DEFINE_double(mapping_keep_ratio, 1.20, "Ratio of the best mapping with the last kept.");
DEFINE_int32(read_pos_max_offset, 10, "Offset tolerance for correct read placement.");

bool Read::read(FILE* fi) {
  static char id[100001];
  static char data[100001];
  static char nesto[100001]; // TODO: sta je treca linija, za sad ignore
  static char kvaliteta[100001]; // TODO: cetvrta valjda kvalitetu oznacava, 
                                 // ignore za sada

  id[100000] = 0;
  if (!fgets(id, sizeof id, fi)) { return 0; }
  assert(id[100000] == 0);

  data[100000] = 0;
  if (!fgets(data, sizeof data, fi)) { return 0; }
  assert(data[100000] == 0);

  nesto[100000] = 0;
  assert(fgets(nesto, sizeof nesto, fi));
  assert(nesto[100000] == 0);

  kvaliteta[100000] = 0;
  assert(fgets(kvaliteta, sizeof kvaliteta, fi));
  assert(kvaliteta[100000] == 0);
  
  this->id_ = id;
  this->data_ = data;
  
  trim(this->id_);
  trim(this->data_);

  return this->id_ != "" && this->data_ != "";
}

void Read::print(FILE* out) {
  fprintf(out, "id: %s\n", id_.c_str());
  fprintf(out, "read: %s\n", data_.c_str());

  fprintf(out, "mappings:\n");
  for (auto mapping : topMappings_) {
    mapping.print(out);
  }
}


int Read::validateWgsimMapping() {
  if (topMappings_.size() == 0) {
    return -2; // nema tog kmera uopce
  }


  vector<string> tokens = Split(id_, '_');
  int pos1 = -1000000, pos2 = -1000000;
  
  string readInGene = tokens[0].substr(1);

  assert(sscanf(tokens[1].c_str(), "%d", &pos1) == 1);
  assert(sscanf(tokens[2].c_str(), "%d", &pos2) == 1);
  
  --pos1; --pos2;

  int i = 0;
  for (auto mapping = topMappings_.rbegin(); 
       mapping != topMappings_.rend(); ++mapping, ++i) {
    bool geneMatch = false;

    // ovo nije nuzno savrseno tocno, ali mislim da se u stvarnosti ne
    // dogadja slucaj kada ne radi
     if (mapping->geneDescriptor.find(readInGene) != string::npos) {
      geneMatch = true;
    }

    if (geneMatch &&
        (abs(mapping->geneBegin-pos1) < FLAGS_read_pos_max_offset ||
         abs(mapping->geneBegin-(pos2-size())) < FLAGS_read_pos_max_offset)) {
      return i;
    }
  }

  // pogresno smjestio
  return -1;
}

void Read::updateMapping(double score, int geneBegin, int geneEnd,
			 int isRC, int geneIdx,
			 string geneDescriptor, 
			 string geneSegment) {
  topMappings_.insert(OneMapping(score, geneBegin, geneEnd, isRC, geneIdx,
				 geneDescriptor, geneSegment));
  
  while (topMappings_.size() > 1 && 
	 topMappings_.rbegin()->score / topMappings_.begin()->score > 
	 FLAGS_mapping_keep_ratio) {
    topMappings_.erase(topMappings_.begin());
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
