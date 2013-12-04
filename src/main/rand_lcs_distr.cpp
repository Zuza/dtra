// For running simulations with various parameter configurations,
// use src/scripts/run_simulations.sh

#include <unistd.h>

#include <gflags/gflags.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "core/klcs.h"
#include "core/ThreadPool.h"
using namespace std;

DEFINE_int32(simulation_threads, sysconf(_SC_NPROCESSORS_ONLN),
	     "Threads used for simulation, defaulting to the number of cores.");
DEFINE_int32(read_len, 1000, "Length of the reads");
DEFINE_int32(simulation_runs, 1000, "Number of performed simulations.");
DEFINE_int32(seed_len, 1, "Seed length");
DEFINE_double(p_err, -1.0, "If set to -1 then rand-to-rand strings are aligned, otherwise rand-to-modified-copy simulations are performed.");
DEFINE_bool(test_all_implementations, false, "");

int min3(int a, int b, int c) { return min(min(a,b),c); }
int max3(int a, int b, int c) { return max(max(a,b),c); }

string create_var_suffix() {
  ostringstream suffix;
  suffix << "_readlen" << FLAGS_read_len 
	 << "_seedlen" << FLAGS_seed_len;
  
  if (FLAGS_p_err < 0) {
    suffix << "_randrand";
  } else {
    suffix << "_perr" << static_cast<int>(FLAGS_p_err*100);
  }
  return suffix.str();
}

const string kNuc = "ACTG";

double p(char a) {
  if (a == 'A') return 0.3;
  if (a == 'C') return 0.2;
  if (a == 'G') return 0.2;
  return 0.3;
}

double p(int i) {
  return p(kNuc[i]);
}

char get_random_base() {
  double r = 1.0*rand()/RAND_MAX;
  for (int j = 0; j < 4; ++j) {
    if (r <= p(j)) {
      return kNuc[j];
    } else {
      r -= p(j);
    }
  }
  assert(0);
}

string generate_string() {
  string ret;
  for (int i = 0; i < FLAGS_read_len; ++i) {
    ret += get_random_base();
  }
  assert(ret.size() == FLAGS_read_len);
  return ret;
}

string generate_similar(const string& a, const double& p_err) {
  string b = a;
  for (int i = 0; i < FLAGS_read_len; ++i) {
    if (1.0*rand()/RAND_MAX <= p_err) {
      b[i] = get_random_base();
    }
  }
  return b;
}

int seeded_lcs(const string& a, const string& b, const int K) {
  int klcs_length = 0;
  if (FLAGS_test_all_implementations) {
    int klcs_length_slow = 0;
    klcs_slow(a, b, K, &klcs_length_slow);
    klcs(a, b, K, &klcs_length);
    // printf("a=%s b=%s\n", a.c_str(), b.c_str());
    // printf("klcs=%d slow=%d\n",
    // 	   klcs_length, klcs_length_slow);
    assert(klcs_length == klcs_length_slow);
  } else {
    klcs(a, b, K, &klcs_length);
  }
  
  return klcs_length;
}

int run_one_simulation() {
  string a = generate_string();
  string b = "";
  
  if (FLAGS_p_err < 0) {
    b = generate_string();
  } else {
    b = generate_similar(a, FLAGS_p_err);
  }

  return seeded_lcs(a, b, FLAGS_seed_len);
}

void output_python_import() {
  puts("from pylab import *"); 
  puts("");
}

void output_python_vector(const string& name, 
			  const vector<double>& vec) {
  printf("%s = [", name.c_str());
  for (int i = 0; i < vec.size(); ++i) {
    if (i) printf(",");
    printf(" %lf", vec[i]);
  }
  puts(" ]");
  puts("");
}

void output_python_plot(const string& plot_name,
			const string& x, const string& y) {
  puts("figure()");
  
  printf("%s,=plot(%s,%s)\n", 
	 plot_name.c_str(), x.c_str(), y.c_str());
  puts("xlabel('LCS Length')");
  puts("ylabel('Probability')");
  //  puts("ylim(-0.1, 1.1)");
  printf("legend([%s], ['%s'], loc=2)\n",
	 plot_name.c_str(), plot_name.c_str());
  puts("show()");
  puts("");
}


int main(int argc, char* argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  //  srand(time(NULL));
  srand(1603);
  
  ThreadPool pool(FLAGS_simulation_threads);
  vector<future<int> > results;
  
  for (int i = 0; i < FLAGS_simulation_runs; ++i) {
    results.push_back(pool.enqueue<int>([] {
	  return run_one_simulation();
	}));
  }
  
  map<int, double> distr;
  for (int i = 0; i < FLAGS_simulation_runs; ++i) {
    distr[results[i].get()] += 1.0/FLAGS_simulation_runs;
  }
  
  vector<double> x;
  vector<double> y;
  double sum_prob = 0;
  
  for (int i = 0; i <= FLAGS_read_len; ++i) {
    double p = distr[i];
    sum_prob += p;
    
    x.push_back(1.0*i/FLAGS_read_len);
    y.push_back(p);
  }

  assert(0.99999 <= sum_prob <= 1.00001);

  string var_suffix = create_var_suffix();
  string plot_name = "plot"+var_suffix;
  string x_var_name = "x"+var_suffix;
  string y_var_name = "y"+var_suffix;
  output_python_import();
  output_python_vector(x_var_name, x);
  output_python_vector(y_var_name, y);
  output_python_plot(plot_name, x_var_name, y_var_name);
  return 0;
}
