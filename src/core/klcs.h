#include <string>
#include <utility>
#include <vector>

// Variant when match pairs are preprocessed, used in the aligner.
void klcs(const std::vector<std::pair<int, int> >& matches,
	  const int k, int* klcs_length,
	  std::vector<std::pair<int, int> >* klcs_reconstruction);

void klcs(const char* aString, size_t aLen, const char* bString, size_t bLen,
	  const int k, int* klcs_length,
	  std::vector<std::pair<int, int> >* klcs_reconstruction);

void klcs_sparse_slow(const char* aString, size_t aLen, const char* bString, size_t bLen,
		      const int k, int* klcs_length,
		      std::vector<std::pair<int, int> >* klcs_reconstruction);

void klcs_sparse_fast(const char* aString, size_t aLen, const char* bString, size_t bLen,
		      const int k, int* klcs_length,
		      std::vector<std::pair<int, int> >* klcs_reconstruction);

bool valid_klcs(const char* aString, size_t aLen, const char* bString, size_t bLen,
		const int k, const int klcs_len,
		const std::vector<std::pair<int, int> >& klcs_recon);

void klcs_slow(const char* aString, size_t aLen, const char* bString, size_t bLen, const int K,
	       int* klcs_length);
