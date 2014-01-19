#include <string>
#include <utility>
#include <vector>

// Variant when match pairs are preprocessed,
// used in the aligner.
void klcs(const std::vector<std::pair<int, int> >& matches,
	  const int k, int* klcs_length,
	  std::vector<std::pair<int, int> >* klcs_reconstruction);

void klcs(const std::string& a, const std::string& b, 
	  const int k, int* klcs_length,
	  std::vector<std::pair<int, int> >* klcs_reconstruction);

void klcs_sparse_slow(const std::string& a, const std::string& b, 
		      const int k, int* klcs_length,
		      std::vector<std::pair<int, int> >* klcs_reconstruction);

void klcs_sparse_fast(const std::string& a, const std::string& b, 
		      const int k, int* klcs_length,
		      std::vector<std::pair<int, int> >* klcs_reconstruction);

bool valid_klcs(const std::string& a, const std::string& b,
		const int k, const int klcs_len,
		const std::vector<std::pair<int, int> >& klcs_recon);

// TODO(fpavetic): Rekonstrukcija
void klcs_slow(const std::string& a, const std::string& b, const int K,
	       int* klcs_length);
