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
