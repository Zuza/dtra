#include <cassert>
#include <cstring>

#include <fstream>
#include <string>
#include <vector>

#include "read.h"
#include "util.h"

using namespace std;

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
  if (fscanf(fi, "%s", data) != 1) return 0;
  assert(data[100000] == 0);

  nesto[100000] = 0;
  if (fscanf(fi, "%s", nesto) != 1) return 0;
  assert(nesto[100000] == 0);

  kvaliteta[100000] = 0;
  if (fscanf(fi, "%s", kvaliteta) != 1) return 0;
  assert(kvaliteta[100000] == 0);
  
  this->id_ = id;
  this->data_ = data;
  
  return 1;
}

void Read::print() {
  printf("id: %s\n", id_.c_str());
  printf("read: %s\n", data_.c_str());

  printf("mappings:\n");
  for (int i = 0; i < topMappings_.size(); ++i) {
    topMappings_[i].print();
  }
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
