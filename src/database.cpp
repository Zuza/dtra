#include "database.h"
using namespace std;

const size_t kMaxBlockSize = 10000000;

Database::Database(const string& databasePath) {
  inputFilePointer_ = fopen(databasePath.c_str(), "rt");
}

Database::~Database() {
  fclose(inputFilePointer_);
}

bool Database::readNextBlock() {
  species_.clear();
  size_t blockSize = 0;

  while (true) {
    Genome g;
    if (!readGenome(&g, inputFilePointer_)) {
      break;
    }

    species_.push_back(g);
    blockSize += g.data_.size();
    if (blockSize > kMaxBlockSize) {
      break;
    }
  }

  return species_.size() > 0;
}
