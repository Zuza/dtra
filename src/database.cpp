#include "database.h"
#include <cstdlib>
using namespace std;

const size_t kMaxBlockSize = 10000000;
const int kSeedLen = 20;
const int kMaxHits = 1000000;

Database::Database(const string& databasePath) {
  inputFilePointer_ = fopen(databasePath.c_str(), "rt");
  if (!inputFilePointer_) {
    fprintf(stderr, "Failed to read database!");
    exit(1);
  }
}

Database::~Database() {
  fclose(inputFilePointer_);
}

bool Database::readNextBlock() {
  species_.clear();
  speciesIndex_.clear();
  currentBlockNoBytes_ = 0;

  while (true) {
    Genome g;
    if (!readGenome(&g, inputFilePointer_)) {
      break;
    }

    species_.push_back(g);
    speciesIndex_.push_back(Index(&g,
    				  kSeedLen,
    				  kMaxHits));

    currentBlockNoBytes_ += g.name_.size();
    currentBlockNoBytes_ += g.data_.size();
    if (currentBlockNoBytes_ > kMaxBlockSize) {
      break;
    }
  }

  return species_.size() > 0;
}

void Database::printNames() {
  for (auto x : species_) {
    printf("%s\n", x.name_.c_str());
  }
}
