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

#include "core/gene.h"
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
DEFINE_string(dna_input_file, "", "Input file from which DNA segments are sampled. If this is not set, DNA segments will be randomly generated.");

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

string generate_string(const int& len) {
  string ret;
  for (int i = 0; i < len; ++i) {
    ret += get_random_base();
  }
  assert(ret.size() == len);
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

struct DnaSampler {
  DnaSampler() {
  }

  void init_from_file(const string& path) {
    FILE* fd = fopen(path.c_str(), "rt");
    assert(fd);

    while (true) {
      shared_ptr<Gene> gene(new Gene());
      if (!readGene(gene.get(), fd)) { break; }
      genes.push_back(gene);
    }
    fclose(fd);
  }

  string remove_ns(const string& a) {
    string ret = "";
    for (int i = 0; i < (int)a.size(); ++i) {
      if (a[i] != 'N') ret += a[i];
      else ret += get_random_base();
    }
    return ret;
  }

  string sample_from_genome(const int& len) {
    long long total_size = 0;
    for (int i = 0; i < genes.size(); ++i) {
      assert(genes[i]->dataSize() >= len);
      total_size += genes[i]->dataSize();
    }
    
    double coin = 1.0*rand()/RAND_MAX; // biramo kromosom prop. s duljinom
    long long running_size = 0;

    for (int i = 0; i < genes.size(); ++i) {
      running_size += genes[i]->dataSize();
      if (1.0*running_size/total_size >= coin) {
	// znam iz kojeg entry-a ulaznog filea uzimam
	int start = rand()%(genes[i]->dataSize()-len+1);
	assert(0 <= start);
	return remove_ns(string(genes[i]->data()+start,
				genes[i]->data()+start+len));
      }
    }

    assert(0);
  }

  pair<string,string> get_random_dna_pair(const double& p_err) {
    pair<string, string> random_pair;

    if (!genes.empty()) { // sampliranje iz genoma
      assert(FLAGS_dna_input_file != "");
      random_pair.first = sample_from_genome(FLAGS_read_len);
      if (p_err < 0) {
	random_pair.second = sample_from_genome(FLAGS_read_len);
      } else {
	random_pair.second = generate_similar(random_pair.first, p_err);
      }
    } else {
      assert(FLAGS_dna_input_file == "");
      random_pair.first = generate_string(FLAGS_read_len);
      if (p_err < 0) {
	random_pair.second = generate_string(FLAGS_read_len);
      } else {
	random_pair.second = generate_similar(random_pair.first, p_err);
      }
    }
    return random_pair;
  }

  vector<shared_ptr<Gene>> genes;
};

DnaSampler global_dna_sampler;

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

int seeded_lcs(const string& a, const string& b, const int K) {
  int klcs_length = 0;

  if (FLAGS_test_all_implementations) {
    klcs_slow(a, b, K, &klcs_length);

    int klcs_sparse_slow_len = 0;
    vector<pair<int, int> > klcs_sparse_slow_recon;
    klcs_sparse_slow(a, b, K, &klcs_sparse_slow_len, &klcs_sparse_slow_recon);

    int klcs_sparse_fast_len = 0;
    vector<pair<int, int> > klcs_sparse_fast_recon;
    klcs_sparse_fast(a, b, K, &klcs_sparse_fast_len, &klcs_sparse_fast_recon);

    printf("%d %d %d\n", klcs_length, klcs_sparse_slow_len,
	   klcs_sparse_fast_len);

    assert(klcs_length == klcs_sparse_slow_len);
    assert(klcs_length == klcs_sparse_fast_len);
    assert(valid_klcs(a, b, K, klcs_sparse_slow_len, 
		      klcs_sparse_slow_recon));
    assert(valid_klcs(a, b, K, klcs_sparse_fast_len,
    		      klcs_sparse_fast_recon));
  } else {
    klcs(a, b, K, &klcs_length, NULL);
  }
  
  return klcs_length;
}

int run_one_simulation() {
  pair<string, string> ab = 
    global_dna_sampler.get_random_dna_pair(FLAGS_p_err);
  const string& a = ab.first;
  const string& b = ab.second;
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

void calculate_distribution(map<int, double>& distr) {
  ThreadPool pool(FLAGS_simulation_threads);
  vector<future<int> > results;
  
  for (int i = 0; i < FLAGS_simulation_runs; ++i) {
    results.push_back(pool.enqueue<int>([] {
	  return run_one_simulation();
	}));
  }
  
  distr.clear();
  for (int i = 0; i < FLAGS_simulation_runs; ++i) {
    distr[results[i].get()] += 1.0/FLAGS_simulation_runs;
  }
}

void fill_xy(map<int, double>& distr, 
	     vector<double>& x, 
	     vector<double>& y) {
  double sum_prob = 0;
  double e_lcs = 0;

  for (int i = 0; i <= FLAGS_read_len; ++i) {
    double p = distr[i];
    sum_prob += p;
    e_lcs += p*i;
    
    x.push_back(1.0*i);
    y.push_back(p);
  }

  assert(0.99999 <= sum_prob <= 1.00001);
  printf("Expected k-LCS=%0.3lf\n", e_lcs);
}
	     

int main(int argc, char* argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  srand(1603);

  if (FLAGS_dna_input_file != "") {
    global_dna_sampler.init_from_file(FLAGS_dna_input_file);
    printf("Initialized dna sampler from file (%s).\n",
	   FLAGS_dna_input_file.c_str());
  }

  map<int, double> distr;
  calculate_distribution(distr);
    
  vector<double> x, y;
  fill_xy(distr, x, y);

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
