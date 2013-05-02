#include "core/index.h"
#include "core/bioinf_util.h"
#include <ctime>
#include <cstdlib>
#include <algorithm>

#include "BufferedBinaryWriter.h"

// for parallel sorting
#include <parallel/algorithm>

using namespace std;

DEFINE_double(avg_multiplier, 5.0, "Index will discard every kmer that " \
              "appears more than avg_freq * avg_multiplier times");
DEFINE_bool(use_parallel_sort, true, "Sorting the index will use " \
            "__gnu_parallel::stable_sort");

Index::iterator::iterator(const Index* idx,
			  const hash_t& hash,
			  const int& querySeedLen) : idx_(idx) {
  // TODO: trenutno je querySeedLen ignoriran, treba
  // dodati da u equal_range usporedbama uzimam samo
  // najznacajnijih querySeedLen znamenki u bazi 4

  // zakomentirani dijelovi micu binary search 
  // kod svakog seeda (jer trazi u kojem se genomu nalazi)
  // zakomentirano je jer nisam spavao i ne da mi se testirati

  Index::Entry tmp; tmp.hash = hash;
  auto pair_lb_ub = 
    equal_range(idx->index_.begin(), idx->index_.end(), tmp);
  begin = distance(idx->index_.begin(), pair_lb_ub.first);
  end   = distance(idx->index_.begin(), pair_lb_ub.second);
  reset();
}

// void Index::iterator::setStartingPos(size_t where) {
//   int position = idx_->index_[where].position;
  
//   auto it = upper_bound(idx_->geneStartingPos_.begin(), 
// 			idx_->geneStartingPos_.end(), position);
//   assert(it != idx_->geneStartingPos_.begin());
//   --it;
//   currStartingPos = distance(idx_->geneStartingPos_.begin(), it);  
// }

void Index::iterator::reset() {
  curr = begin;
  // setStartingPos(curr);
}

bool Index::iterator::done() {
  return curr == end;
}

void Index::iterator::advance() {
  ++curr;

  // size_t next = curr+1;
  // setStartingPos(curr = next);
  
  // if (next < end) {
  //   size_t currPosition = idx_->index_[curr].position;
  //   size_t nextPosition = idx_->index_[next].position;

  //   if (nextPosition < currPosition) {
  //     // ovo je potrebno u slucaju kad je querySeedLen manji
  //     // od duljine kojom smo izgradili index. tad se moze
  //     // dogoditi da unutar istog hash ne budu sortirane pozicije
  //     // vec se sastoje od slijepljenih sortiranih nizova
  //     setStartingPos(next); 
  //   } else {
  //     if (currStartingPos+1 < idx_->geneStartingPos_.size()) {
  // 	if (nextPosition >= idx_->geneStartingPos_[currStartingPos+1]) {
  // 	  ++currStartingPos;
  // 	}
  //     }
  //   }
  // }

  // curr = next;
}

pair<unsigned int, unsigned int> Index::iterator::get() {
  return idx_->position_to_gene_position(idx_->index_[curr].position);
  // int position = 
  //   idx_->index_[curr].position - 
  //   idx_->geneStartingPos_[currStartingPos];
  // return make_pair<unsigned int, unsigned int>(currStartingPos,
  // 					       position);
}

Index::Index(int seedLength) : seedLength_(seedLength) {
  indexPrepared_ = false;
  startingPos_ = 0;
  assert(seedLength_ <= 32); 
}

void Index::insertGene(Gene* gene) { 
  // TODO: prealociraj sve!

  hash_t subtractPower = 1;
  for (int i = 0; i < seedLength_; ++i) {
    subtractPower *= 4;
  }

  hash_t rollingHash = 0;

  for (size_t i = 0; i < gene->dataSize(); ++i) {
    int next = baseToInt(gene->data()[i]);

    rollingHash = rollingHash * 4 + next;
    if (i >= seedLength_) {
      rollingHash %= subtractPower; // TODO(zuza): ANDanje umjesto moda
    }

    if (i+1 >= seedLength_) {
      Index::Entry tmp;
      tmp.position = startingPos_ + i+1 - seedLength_;
      tmp.hash = rollingHash;
      index_.push_back(tmp);
    }
  }

  geneStartingPos_.push_back(startingPos_);
  startingPos_ += gene->dataSize();
}

void Index::discardFrequentSeeds() {
  assert(indexPrepared_);

  hash_t last_hash = index_[0].hash + 1;
  long long run = 0;

  long long run_max = 0;
  long long run_sum = 0;
  long long run_cnt = 0;
  long long kmers_cnt = 0;

  for (auto entry : index_) {
    ++kmers_cnt;

    if (last_hash == entry.hash) {
      ++run;
    } else if (run > 0) {
      /////
      if (run > run_max) run_max = run;
      run_sum += run;
      ++run_cnt;
      /////
      run = 1;
      last_hash = entry.hash;
    }
  }
  if (run > 0) {
    /////
    if (run > run_max) run_max = run;
    run_sum += run;
    ++run_cnt;
    /////
  }

  double avg_run = double(run_sum) / run_cnt;
  printf("Maximum k-mer frequency = %lld\n", run_max);
  printf("Average k-mer frequency = %.2lf (discarding freqs larger than %.2lf)\n", avg_run, FLAGS_avg_multiplier * avg_run);

  long long runs_left = 0;
  long long kmers_left = 0;

  size_t write_over = 0;
  for (size_t i = 0; i < index_.size(); ) {
    hash_t last_hash = index_[i].hash;
    long long run = 0;

    size_t j;
    for (j = i; j < index_.size(); ++j) {
      if (last_hash == index_[j].hash)
        ++run;
      else
        break;
    }

    // if (run <= FLAGS_avg_multiplier * run_avg)
    if (run_cnt * run <= FLAGS_avg_multiplier * run_sum) {
      ++runs_left;
      kmers_left += run;

      for (long long it = 0; it < run; ++it) {
        index_[write_over] = index_[i + it];
        ++write_over;
      }
    }

    i = j;
  }
  index_.resize(write_over);

  printf("Different hashes left = %lld (%.2lf%%)\n", runs_left, double(runs_left) / run_cnt * 100.0);
  printf("Kmers left = %lld (%.2lf%%)\n", kmers_left, double(kmers_left) / kmers_cnt * 100.0);
}


pair<unsigned int, unsigned int> Index::position_to_gene_position(size_t position) const {
  auto it = upper_bound(geneStartingPos_.begin(), geneStartingPos_.end(), position);
  assert(it != geneStartingPos_.begin());
  --it;
  return make_pair<unsigned int, unsigned int>(it - geneStartingPos_.begin(), position - *it);
}

Index::iterator Index::getPositions(const hash_t& hash,
				    const int& querySeedLen) {
#ifdef DEBUG
  assert(indexPrepared_);
#endif
  return iterator(this, hash, querySeedLen);
}

void Index::prepareIndex() {
  // positions are in increasing order
  if (FLAGS_use_parallel_sort) {
    __gnu_parallel::stable_sort(index_.begin(), index_.end()); // STABLE!!
  } else {
    std::stable_sort(index_.begin(), index_.end());
  }
  indexPrepared_ = true;
}

void Index::writeIndex(BufferedBinaryWriter& writer) {
  assert(indexPrepared_);
  writer.writeSigned32(seedLength_);
  writer.writeVector(geneStartingPos_);
  writer.writeVector(index_);
}

void Index::readIndex(BufferedBinaryReader& reader) {
  indexPrepared_ = true;
  int storedSeedLength; reader.readSigned32(&storedSeedLength);
  assert(seedLength_ == storedSeedLength);
  reader.readVector(geneStartingPos_);
  reader.readVector(index_);
}
