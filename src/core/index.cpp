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

pair<unsigned int, unsigned int> Index::position_to_gene_position(size_t position) {
  auto it = upper_bound(geneStartingPos_.begin(), geneStartingPos_.end(), position);
  assert(it != geneStartingPos_.begin());
  --it;
  return make_pair<unsigned int, unsigned int>(it - geneStartingPos_.begin(), position - *it);
}

void Index::getPositions(vector<pair<unsigned int, unsigned int> >* retVal, hash_t hash) {
  assert(indexPrepared_);
  assert(retVal->empty());
  Entry tmp; tmp.hash = hash;
  auto pair_lb_ub = equal_range(index_.begin(), index_.end(), tmp);
  for ( ; pair_lb_ub.first != pair_lb_ub.second; ++pair_lb_ub.first) {
    Entry& entry = *pair_lb_ub.first;
    retVal->push_back(position_to_gene_position(entry.position));
  }
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
