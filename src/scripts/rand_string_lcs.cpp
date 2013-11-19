#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>
#include <vector>
using namespace std;

int max3(int a, int b, int c) { return max(max(a,b),c); }

const int N = 1000;
const int ITERACIJA = 1000;
map<int, double> distr;
int dp[N+1][N+1];

string nuc = "ACTG";
double p(char a) {
  if (a == 'A') return 0.3;
  if (a == 'C') return 0.2;
  if (a == 'G') return 0.2;
  return 0.3;
}

double p(int i) {
  return p(nuc[i]);
}

char get_random_base() {
  double r = 1.0*rand()/RAND_MAX;
  for (int j = 0; j < 4; ++j) {
    if (r <= p(j)) {
      return nuc[j];
    } else {
      r -= p(j);
    }
  }
  assert(0);
}

string generate_string() {
  string ret;
  for (int i = 0; i < N; ++i) {
    ret += get_random_base();
  }
  assert(ret.size() == N);
  return ret;
}

string generate_similar(const string& a, const double& p_err) {
  string b = a;
  for (int i = 0; i < N; ++i) {
    if (1.0*rand()/RAND_MAX <= p_err) {
      b[i] = get_random_base();
    }
  }
  return b;
}

int lcsK(const string& a, const string& b, const int K) {
  for (int i = 0; i <= N; ++i) dp[i][0] = 0;
  for (int j = 0; j <= N; ++j) dp[0][j] = 0;

  for (int i = 1; i <= N; ++i) {
    for (int j = 1; j <= N; ++j) {
      dp[i][j] = max(dp[i-1][j], dp[i][j-1]);

      int lo = K, hi = min(i,j);

      while (lo <= hi) {
	int k = (lo+hi)/2;

	if (a.substr(i-k, k) == b.substr(j-k, k)) {
	  lo = k+1;
	  dp[i][j] = max(dp[i][j], dp[i-k][j-k]+k);
	} else {
	  hi = k-1;
	}
      }
    }
  }

  return dp[N][N];
}

int f() {
  string a = generate_string();
  string b = generate_string();
  //string b = generate_similar(a, 0.10);
  return lcsK(a,b,1);
}

void output_matlab_vector(const string& name, const vector<double>& vec) {
  printf("%s = [", name.c_str());
  for (int i = 0; i < vec.size(); ++i) {
    printf(" %lf", vec[i]);
  }
  puts(" ];");
}

int main(void) {
  for (int i = 0; i < ITERACIJA; ++i) {
    distr[f()] += 1.0/ITERACIJA;
  }

  vector<double> x;
  vector<double> y;
  double sum_prob = 0;

  for (int i = 0; i <= N; ++i) {
    double p = distr[i];
    printf("%5d: %0.10lf\n", i, p);
    sum_prob += p;

    x.push_back(i);
    y.push_back(p);
  }
  
  printf("Sum of probabilities = %0.2lf\n", sum_prob);

  output_matlab_vector("x", x);
  output_matlab_vector("y", y);
  return 0;
}
