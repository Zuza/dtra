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


static void klcs_base(const int k, const vector<pair<int, int> >& matches,
		      int* klcs_length) {
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
}

// Tocna, ali slozenija implementacija.
// 
// static void klcs_base(const int k, const vector<pair<int, int> >& matches,
// 		      int* klcs_length) {
//   struct PointOfInterest {
//     int i, j;
//     enum Type {BEGIN = 0, END = 1};
//     Type type;

//     PointOfInterest(int i, int j, Type type) : 
//       i(i), j(j), type(type) {}
    
//     bool operator < (const PointOfInterest& other) const {
//       if (i != other.i) { return i < other.i; }
//       if (j != other.j) { return j < other.j; }
//       return type < other.type;
//     }

//     bool operator == (const PointOfInterest& other) const {
//       return 
// 	i == other.i && 
// 	j == other.j && 
// 	type == other.type;
//     }
//   };

//   struct DiagonalDpVal {
//     int secondary_diagonal;
//     int value;

//     DiagonalDpVal() {}

//     DiagonalDpVal(int secondary_diagonal, int value):
//       secondary_diagonal(secondary_diagonal), value(value) {}

//     bool operator < (const DiagonalDpVal& other) const {
//       return value < other.value;
//     }

//     bool operator == (const DiagonalDpVal& other) const {
//       return 
// 	secondary_diagonal == other.secondary_diagonal &&
// 	value == other.value;
//     }
//   };


//   vector<PointOfInterest> poi;
//   for (int i = 0; i < matches.size(); ++i) {
//     poi.push_back(PointOfInterest(matches[i].first,
// 				  matches[i].second,
// 				  PointOfInterest::Type::BEGIN));
//     poi.push_back(PointOfInterest(matches[i].first+k-1,
// 				  matches[i].second+k-1,
// 				  PointOfInterest::Type::END));
//   }

//   sort(poi.begin(), poi.end());
//   // ovaj unique mozda nije potreban?
//   poi.erase(unique(poi.begin(), poi.end()), poi.end());

//   int n = 0;
//   for (auto it = poi.begin(); it != poi.end(); ++it) {
//     n = max(n, it->i+1);
//     n = max(n, it->j+1);
//   }

//   *klcs_length = 0;
//   // Indexed by column:
//   FenwickMax<int> first_columns_max(n);

//   vector<pair<int, int> > curr_row_vals;
//   int curr_row = -1;
//   int curr_row_max_val = 0;

//   vector<MonotonicQueue<DiagonalDpVal> > diag_lcs_k(2*n);
//   vector<MonotonicQueue<DiagonalDpVal> > diag_lcs_2k(2*n);

//   for (auto it = poi.begin(); it != poi.end(); ++it) {
//     int i = it->i;
//     int j = it->j;
//     PointOfInterest::Type type = it->type;

//     if (i != curr_row) {
//       while (!curr_row_vals.empty()) {
//   	int col = curr_row_vals.back().first;
//   	int val = curr_row_vals.back().second;
//   	curr_row_vals.pop_back();
//   	first_columns_max.update(col, val);
//       }
//       curr_row = i;
//       curr_row_max_val = 0;
//     }

//     int primary_diagonal = n-1+i-j;
//     int secondary_diagonal = i+j;
    
//     MonotonicQueue<DiagonalDpVal>& lcs_k =
//       diag_lcs_k[primary_diagonal];
//     MonotonicQueue<DiagonalDpVal>& lcs_2k =
//       diag_lcs_2k[primary_diagonal];

//     while (!lcs_k.empty() && 
// 	   (secondary_diagonal-lcs_k.front().secondary_diagonal)/2 >= 
// 	   k - (PointOfInterest::Type::END == type)) {
//       lcs_2k.push(lcs_k.front());
//       lcs_k.pop();
//     }
//     while (!lcs_2k.empty() && 
// 	   (secondary_diagonal-lcs_2k.front().secondary_diagonal)/2 >= 2*k) {
//       lcs_2k.pop();
//     }
    
//     if (type == PointOfInterest::Type::BEGIN) {
//       // 1) nastavak na neki zatvoreni
//       int lcs_ij = (j>0?first_columns_max.get(j-1):0);
    
//       //printf("IZ LOGARITAMSKE=%d\n", lcs_ij);
//       if (!lcs_k.empty()) {
// 	int prev_val = lcs_k.max().value;
// 	//printf("BEGIN prev_val=%d\n", prev_val);
// 	lcs_ij = max(lcs_ij, i - 1 + prev_val);
// 	//printf("UPDATED lcs_ij=%d\n", lcs_ij);
//       }

//       // printf("PUSHAM u lcs_k na diag=%d vrijednost=%d\n", 
//       // 	     primary_diagonal, lcs_ij - (i-1));
//       lcs_k.push(DiagonalDpVal(secondary_diagonal, lcs_ij - (i - 1)));
//       // printf("%s: %d %d = %d\n", 
//       // 	     type == PointOfInterest::Type::BEGIN ? "BEGIN": "END",
//       // 	     i, j, lcs_ij - (i-1));
//       // fflush(stdout);
//     } else {
//       int lcs_ij = 0;
//       // printf("Nalazim se na %d-toj glavnoj dijagonali.\n",
//       // 	     primary_diagonal);
//       // printf("lcs_k.size()=%d, lcs_2k.size()=%d\n",
//       // 	     lcs_k.size(), lcs_2k.size());

//       if (!lcs_2k.empty()) {
// 	int prev_val = lcs_2k.max().value;
// 	// printf("END prev_val=%d\n", prev_val);
// 	lcs_ij = max(lcs_ij, i + prev_val);
//       }
      
//       curr_row_vals.push_back(make_pair(j, lcs_ij));
//       curr_row_max_val = max(curr_row_max_val, lcs_ij);
//       *klcs_length = max(*klcs_length, lcs_ij);

//       // printf("%s: %d %d = %d\n", 
//       // 	     type == PointOfInterest::Type::BEGIN ? "BEGIN": "END",
//       // 	     i, j, lcs_ij);
//       // fflush(stdout);
//     }
//   }
// }

void klcs(const string& a, const string& b,
	  const int k, int* klcs_length) {
  vector<pair<int, int> > matches;
  get_matches(a, b, k, &matches);
  klcs_base(k, matches, klcs_length);
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
