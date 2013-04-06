#include "database.h"
#include <cstdlib>
using namespace std;

const size_t kMaxBlockSize = 10000000;
const int kSeedLen = 20;

Database::Database(const string& databasePath,
		   const string& indexFilePath,
		   const bool createIndex) : createIndex_(createIndex) {
  dbFilePointer_ = fopen(databasePath.c_str(), "rt");
  if (!dbFilePointer_) {
    fprintf(stderr, "Failed to read database!\n");
    exit(1);
  }
  
  // TODO: mozda kasnije vidjeti je li potrebna ideja ili
  //       korisno imati index zapisan kao binarni file
  indexFilePointer_ = fopen(indexFilePath.c_str(), 
			    createIndex_ ? "wb": "rb");
  if (!indexFilePointer_) {
    fprintf(stderr, "Failed to open index file!\n");
    exit(1);
  }

  if (!createIndex) { // initialize BufferedReader
    bufferedReader_ = shared_ptr<BufferedBinaryReader>(
                          new BufferedBinaryReader(indexFilePointer_));
  } else {
    bufferedWriter_ = shared_ptr<BufferedBinaryWriter>(
		          new BufferedBinaryWriter(indexFilePointer_));
  }
}

Database::~Database() {
  fclose(dbFilePointer_);
  fclose(indexFilePointer_);
}

bool Database::readNextBlock() { 
  species_.clear();
  speciesIndex_.clear();
  currentBlockNoBytes_ = 0;

  for (int iter = 0; ; ++iter) {
    shared_ptr<Genome> g(new Genome());
    if (!readGenome(g.get(), dbFilePointer_)) {
      break;
    }
    species_.push_back(g);

    shared_ptr<Index> in(new Index(kSeedLen));
    if (createIndex_) {
      in->create(g.get());
      in->appendToBufferedBinaryWriter(*bufferedWriter_);
    } else {
      in->readNextFromBufferedReader(*bufferedReader_);
    }
    //printf("%llu\n", in->checksum());
    speciesIndex_.push_back(in);
    currentBlockNoBytes_ += g->name_.size();
    currentBlockNoBytes_ += g->data_.size();
    if (currentBlockNoBytes_ > kMaxBlockSize) {
      break;
    }
  }

  return species_.size() > 0;
}

void Database::printNames() {
  for (auto x : species_) {
    printf("%s\n", x->name_.c_str());
  }
}
