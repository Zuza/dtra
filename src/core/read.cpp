#include <cstring>

#include <fstream>
#include <string>
#include <vector>

#include "read.h"
#include "util.h"

using namespace std;

int Read::readFromStdin() {
  static char id[100000];
  static char data[100000];
  static char nesto[100000]; // TODO: sta je treca linija, za sad ignore
  static char kvaliteta[100000]; // TODO: cetvrta valjda kvalitetu oznacava, 
                                 // ignore za sada

  if (scanf("%s", id) != 1) return 0;
  if (scanf("%s", data) != 1) return 0;
  if (scanf("%s", nesto) != 1) return 0;
  if (scanf("%s", kvaliteta) != 1) return 0;
  
  this->id_ = id;
  this->data_ = data;
  
  return 1;
}

void Read::print() {
  printf("id: %s\n", id_.c_str());
  printf("read: %s\n", data_.c_str());
}
