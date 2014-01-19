#include "bioinf_util.h"
#include "fenwick.h"
#include "klcs.h"
#include "monotonic_queue.h"

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cassert>
using namespace std;

/**
 * <matches> is filled with pairs (i,j) meaning that for 
 * every such pair a[i...i+k-1] == b[j...j+k-1]
 */
static void get_matches(const string& a, const string& b,
			const int k, vector<pair<int, int> >* matches) {
  assert(k < 32);

  if (matches != NULL) {
    matches->clear();
  } else {
    return;
  }

  typedef unordered_multimap<uint64_t, int> MatchIndexType;
  unique_ptr<MatchIndexType> match_index =
    unique_ptr<MatchIndexType>(new MatchIndexType());

  uint64_t hash_mod = (1LL<<(2*k))-1;
  uint64_t rolling_hash = 0;
  for (int i = 0; i < a.size(); ++i) {
    assert(isBase(a[i]));

    rolling_hash = rolling_hash * 4 + baseToInt(a[i]);
    rolling_hash &= hash_mod;

    if (i+1 >= k) {
      match_index->insert(MatchIndexType::value_type(rolling_hash, i-k+1));
    }
  }

  rolling_hash = 0;
  for (int i = 0; i < b.size(); ++i) {
    assert(isBase(b[i]));
    
    rolling_hash = rolling_hash * 4 + baseToInt(b[i]);
    rolling_hash &= hash_mod;

    if (i+1 >= k) {
      auto positions_in_a = match_index->equal_range(rolling_hash);
      for (auto it = positions_in_a.first; it != positions_in_a.second; ++it) {
	matches->push_back(make_pair(it->second, i-k+1));
      }
    }
  }
}

static void fill_klcs_reconstruction(const vector<pair<int, int> >& matches,
				     const int k,
				     const vector<int>& prev_idx,
				     const int last_idx,
				     vector<pair<int, int> >* klcs_recon) {
  klcs_recon->clear();
  for (int i = last_idx; i != -1; i = prev_idx[i]) {
    int r = matches[i].first+k-1;
    int c = matches[i].second+k-1;

    if (prev_idx[i] == -1 || 
	(matches[prev_idx[i]].first + k <= matches[i].first &&
	 matches[prev_idx[i]].second + k <= matches[i].second)) {
      // Uzimam sigurno cijelog
      for (int j = 0; j < k; ++j, --r, --c) {
	klcs_recon->push_back(make_pair(r, c));
      }
    } else { 
      // inace je sigurno nastavak iste primarne dijagonale,
      // uzimam samo nastavak
      int curr = i;
      int prev = prev_idx[i];

      int curr_secondary_diag = (matches[curr].first + matches[curr].second) / 2;
      int prev_secondary_diag = (matches[prev].first + matches[prev].second) / 2;
      assert(prev_secondary_diag < curr_secondary_diag);

      for (int j = prev_secondary_diag; j < curr_secondary_diag; ++j, --r, --c) {
	klcs_recon->push_back(make_pair(r, c));
      }
    }
  }

  reverse(klcs_recon->begin(), klcs_recon->end());
}

static void klcs_sparse_fast(const vector<pair<int, int> >& matches,
			     const int k, int* klcs_length,
			     vector<pair<int, int> >* klcs_reconstruction) {
  typedef pair<int, int> PointOfInterest;

  struct DiagonalDpVal {
    int secondary_diagonal;
    int value;
    
    DiagonalDpVal() {}
    
    DiagonalDpVal(int secondary_diagonal, int value):
      secondary_diagonal(secondary_diagonal), value(value) {}

    bool operator < (const DiagonalDpVal& other) const {
      return value < other.value;
    }

    bool operator == (const DiagonalDpVal& other) const {
      return 
	secondary_diagonal == other.secondary_diagonal &&
	value == other.value;
    }
  };

  vector<PointOfInterest> poi;
  for (int i = 0; i < matches.size(); ++i) {
    poi.push_back(PointOfInterest(matches[i].first+k-1,
				  matches[i].second+k-1));
  }

  sort(poi.begin(), poi.end());

  int n = 0;
  for (auto it = poi.begin(); it != poi.end(); ++it) {
    n = max(n, it->first+1);
    n = max(n, it->second+1);
  }

  *klcs_length = 0;
  // Indexed by column:
  FenwickMax<int> first_columns_max(n);
  queue<std::tuple<int, int, int> > update_queue;
  vector<MonotonicQueue<DiagonalDpVal> > diag_lcs_k(2*n);

  for (auto it = poi.begin(); it != poi.end(); ++it) {
    int i = it->first;
    int j = it->second;

    while (!update_queue.empty() &&
	   get<0>(update_queue.front())+k <= i) {
      auto t = update_queue.front();
      update_queue.pop();
      int col = get<1>(t);
      int val = get<2>(t);
      first_columns_max.update(col, val);
    }

    int primary_diagonal = n-1+i-j;
    int secondary_diagonal = i+j;
    MonotonicQueue<DiagonalDpVal>& lcs_k =
      diag_lcs_k[primary_diagonal];

    while (!lcs_k.empty() && 
	   (secondary_diagonal-lcs_k.front().secondary_diagonal)/2 >= k) {
      lcs_k.pop();
    }
    
    int lcs_ij = (j>=k?first_columns_max.get(j-k):0)+k;
    if (!lcs_k.empty()) {
      int prev_val = lcs_k.max().value;
      lcs_ij = max(lcs_ij, i + prev_val);     
    }

    lcs_k.push(DiagonalDpVal(secondary_diagonal, lcs_ij - i));
    update_queue.push(make_tuple(i,j,lcs_ij));
    *klcs_length = max(*klcs_length, lcs_ij);
  }

  if (klcs_reconstruction) {
    // TODO: rekonstrukcija
  }
}

static void klcs_sparse_slow(vector<pair<int, int> > matches,
			     const int k, int* klcs_length,
			     vector<pair<int, int> >* klcs_reconstruction) {
  sort(matches.begin(), matches.end());
  int n = matches.size();

  if (n == 0) {
    *klcs_length = 0;
    klcs_reconstruction->clear(); 
  } else {
    vector<int> dp(n);
    vector<int> recon(n);
    int best_idx = 0;
    
    for (int i = 0; i < n; ++i) {
      dp[i] = k;
      recon[i] = -1;

      int end_row_i = matches[i].first+k-1;
      int end_col_i = matches[i].second+k-1;

      int primary_diagonal_i = end_row_i - end_col_i;
      int secondary_diagonal_i = (end_row_i + end_col_i)/2;

      for (int j = i-1; j >= 0; --j) {
	// 1) Uzimam cijeli match interval i nastavljam neki
	// match koji je ranije vec zavrsio.
	if (matches[j].first + k <= matches[i].first &&
	    matches[j].second + k <= matches[i].second) {
	  if (dp[j] + k > dp[i]) {
	    dp[i] = dp[j] + k;
	    recon[i] = j;
	  }
	}

	// 2) Nastavak po istoj dijagonali.
	int end_row_j = matches[j].first+k-1;
	int end_col_j = matches[j].second+k-1;
	int primary_diagonal_j = end_row_j - end_col_j;
	int secondary_diagonal_j = (end_row_j + end_col_j)/2;

	if (primary_diagonal_i == primary_diagonal_j &&
	    secondary_diagonal_i > secondary_diagonal_j &&
	    secondary_diagonal_i - secondary_diagonal_j < k) {
	  if (dp[j] + secondary_diagonal_i - secondary_diagonal_j > dp[i]) {
	    dp[i] = dp[j] + secondary_diagonal_i - secondary_diagonal_j;
	    recon[i] = j;
	  }
	}
      }

      if (dp[i] > dp[best_idx]) {
	best_idx = i;
	*klcs_length = dp[i];
      }
    }

    if (klcs_reconstruction) {
      fill_klcs_reconstruction(matches, k, recon, best_idx, 
			       klcs_reconstruction);
    }
  }
}

void klcs_sparse_slow(const string& a, const string& b,
		      const int k, int* klcs_length,
		      vector<pair<int, int> >* klcs_reconstruction) {
  vector<pair<int, int> > matches;
  get_matches(a, b, k, &matches);
  klcs_sparse_slow(matches, k, klcs_length, klcs_reconstruction);  
}

void klcs_sparse_fast(const string& a, const string& b,
		      const int k, int* klcs_length,
		      vector<pair<int, int> >* klcs_reconstruction) {
  vector<pair<int, int> > matches;
  get_matches(a, b, k, &matches);
  klcs_sparse_fast(matches, k, klcs_length, klcs_reconstruction);  
}

void klcs(const vector<pair<int, int> >& matches,
	  const int k, int* klcs_length,
	  vector<pair<int, int> >* klcs_reconstruction) {
  if (matches.size() < 500) {
    klcs_sparse_slow(matches, k, klcs_length, klcs_reconstruction);
  } else {
    klcs_sparse_fast(matches, k, klcs_length, klcs_reconstruction);
  }
}

void klcs(const string& a, const string& b,
	  const int k, int* klcs_length,
	  vector<pair<int, int> >* klcs_reconstruction) {
  vector<pair<int, int> > matches;
  get_matches(a, b, k, &matches);
  klcs(matches, k, klcs_length, klcs_reconstruction);  
}

bool valid_klcs(const string& a, const string& b,
		const int k, const int klcs_len,
		const vector<pair<int, int> >& klcs_recon) {
  // 1) Ensure correct length.
  if (klcs_len != klcs_recon.size()) {
    return false;
  }

  // 2) Ensure chars corresponding to the indices match.
  for (auto match: klcs_recon) {
    int i = match.first;
    int j = match.second;
    
    if (i < 0 || i >= a.size()) { return false; }
    if (j < 0 || j >= b.size()) { return false; }
    if (a[i] != b[j]) { return false; }
  }

  // 3) Ensure runs of indices have at least length of k.
  int run_a = 1;
  int run_b = 1;
  for (size_t i = 1; i < klcs_recon.size(); ++i) {
    if (klcs_recon[i-1].first >= klcs_recon[i].first) { return false; }
    if (klcs_recon[i-1].second >= klcs_recon[i].second) { return false; }

    if (klcs_recon[i-1].first+1 == klcs_recon[i].first) { ++run_a; }
    if (klcs_recon[i-1].second+1 == klcs_recon[i].second) { ++run_b; }

    if (i+1 == klcs_recon.size() || 
	klcs_recon[i-1].first+1 != klcs_recon[i].first) {
      if (run_a < k) { return false; }
      run_a = 1;
    }

    if (i+1 == klcs_recon.size() ||
	klcs_recon[i-1].second+1 != klcs_recon[i].second) {
      if (run_b < k) { return false; }
      run_b = 1;
    }
  }

  return true;
}

static int min3(int a, int b, int c) { return min(min(a,b),c); }

void klcs_slow(const string& a, const string& b, const int K,
	       int* klcs_length) {
  // Svaka dretva ce ponovno alocirati ovo polje. 
  // Tako je napravljeno zbog jednostavnosti. Ovaj proces 
  // cemo ionako pokretati samo malen broj puta,
  // tek toliko da generiramo simulacijske rezultate.
  vector<vector<int> > dp(a.size()+1, 
			  vector<int>(b.size()+1));

  for (int i = 0; i <= a.size(); ++i) dp[i][0] = 0;
  for (int j = 0; j <= b.size(); ++j) dp[0][j] = 0;

  for (int i = 1; i <= a.size(); ++i) {
    for (int j = 1; j <= b.size(); ++j) {
      dp[i][j] = max(dp[i-1][j], dp[i][j-1]);

      // 2*K je mislim dovoljna granica jer ce sve iznad toga
      // biti pokriveno nekim prethodnim seedom
      for (int k = 1; k <= min3(i, j, 2*K); ++k) {
	char aa = a[i-k];
	char bb = b[j-k];
	if (aa != bb) {
	  break;
	}

	if (k >= K) {
	  dp[i][j] = max(dp[i][j], dp[i-k][j-k]+k);
	  //printf("i=%d j=%d k=%d dp=%d\n", i, j, k, dp[i][j]);
	}
      }
    }
  }

  *klcs_length = dp[a.size()][b.size()];
}
