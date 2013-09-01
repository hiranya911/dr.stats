#include <iostream>
#include <fstream>
#include <cstdlib>
#include <limits>

#include "boost/program_options.hpp"
#include "dvect.h"
#include "stats.h"

using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::istream;

void compute_simple_mode(istream & in);
void compute_vector_centroid(istream & in);
void compute_vector_kmeans(istream & in, const int k, const int rounds, const bool verbose);
void compute_vector_xmeans(istream & in, const int k_min, const int k_max, const int rounds, const bool verbose);

int main(int argc, char** argv) {
  namespace po = boost::program_options;

  string mode;
  string input_file;
  int precision = 10;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "Print help message")
    ("mode,m", po::value<string>(&mode)->required(), "Calculation mode")
    ("precision,p", po::value<int>(&precision), "Decimal point precision (for printing on the console)")
    ("input,i", po::value<string>(&input_file), "Input file")
    ("clusters,k", po::value<string>(), "Number of clusters or cluster range (used for kmeans)")
    ("rounds,r", po::value<int>(), "Number of rounds to run the calculation")
    ("verbose,v", "Enable verbose output (only supported by some modes)")
    ;

  po::positional_options_description pos_desc; 
  pos_desc.add("input", 1); 

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
    
    if (vm.count("help")) {
      cout << "Usage: drstats [-h][-p PRECISION][-i INPUT_FILE] -m MODE input_file" << endl;
      cout << desc << endl;
      return 0;
    }

    po::notify(vm);

    if (vm.count("mode")) {
      mode = vm["mode"].as<string>();
    }

    cout.precision(precision);
    istream* in = &cin;
    if (input_file != "") {
      in = new std::ifstream(input_file.c_str(), std::ifstream::in);
      std::ifstream* fin = (std::ifstream*) in;
      if (!(*fin)) {
	cerr << "Failed to open input file: " << input_file << endl;
	delete in;
	return 1;
      }
    }
    bool error = false;

    if (mode == "simple") {
      compute_simple_mode(*in);
    } else if (mode == "cent") {
      compute_vector_centroid(*in);
    } else if (mode == "kmeans") {
      if (vm.count("clusters")) {
	string k_str = vm["clusters"].as<string>();
	int k = atoi(k_str.c_str());
	int rounds = 20;
	if (vm.count("rounds")) {
	  rounds = vm["rounds"].as<int>();
	}
	bool verbose = false;
	if (vm.count("verbose")) {
	  verbose = true;
	}
	compute_vector_kmeans(*in, k, rounds, verbose);
      } else {
	cerr << "Error: k is required but missing\n";
	error = true;
      }
    } else if (mode == "xmeans") {
      if (vm.count("clusters")) {
	string k_str = vm["clusters"].as<string>();
	std::stringstream ss(k_str);
	string item;
	std::vector<string> elems;
	while (std::getline(ss, item, ':')) {
	  elems.push_back(item);
	}
	if ((int) elems.size() != 2) {
	  cerr << "k value must be specified as a range of the form k_min:k_max\n";
	  error = true;
	}

	int k_min = atoi(elems.at(0).c_str());
	int k_max = atoi(elems.at(1).c_str());
	
	int rounds = 20;
	if (vm.count("rounds")) {
	  rounds = vm["rounds"].as<int>();
	}
	bool verbose = false;
	if (vm.count("verbose")) {
	  verbose = true;
	}
	compute_vector_xmeans(*in, k_min, k_max, rounds, verbose);
      } else {
	cerr << "Error: k (range specifier) is required but missing\n";
	error = true;
      }
    } else {
      cerr << "Unsupported calculation mode: " << mode << endl;
      error = true;
    }

    if (input_file != "") {
      std::ifstream * fin = (std::ifstream*) in;
      fin->close(); 
      delete in;
    }
    if (error) {
      return 1;
    }

  } catch (po::required_option & e) {
    if (e.get_option_name() != "--input") {
      cerr << "Error: " << e.what() << endl;
    } else {
      cerr << "Error: input file is required but missing\n";
    }
    cerr << "Exiting...\n";
    return 1;
  } catch (std::exception & e) {
    cerr << "Error: " << e.what() << endl;
    cerr << "Exiting...\n";
    return 1;
  }

  return 0;
}

void compute_simple_mode(istream & in) {
  dvect numbers;
  long size = dvect_load(in, numbers);

  if (size <=  0) {
    cerr << "Failed to load any input data\n";
    return;
  }

  cout << "Number of data points: " << numbers.size() << endl;
  cout << "Total: " << stats_sum(numbers) << endl;
  cout << "Mean (Average): " << stats_mean(numbers) << endl;
  cout << "Std. Deviation: " << stats_stddev(numbers) << endl;
  cout << "Median: " << stats_median(numbers) << endl;

  double min = stats_min(numbers);
  double max = stats_max(numbers);
  cout << "Range: " << max - min << endl;
  cout << "Min: " << min << endl;
  cout << "Max: " << max << endl;
}

void compute_vector_centroid(istream & in) {
  dvectlist vectors;
  long size = dvect_load(in, vectors);

  if (size <= 0) {
    cerr << "Failed to load any input data\n";
    return;
  }

  dvect result = stats_vector_centroid(vectors);
  cout << "Number of vectors: " << size << endl;
  cout << "Vector dimensions: " << vectors.front().size() << endl;
  cout << "Centroid: ";
  dvect_print(result);
}

void compute_vector_kmeans(istream & in, const int k, const int rounds, const bool verbose) {
  if (k < 1) {
    cerr << "Cluster count (k) must be greater than 0\n";
    return;
  } else if (rounds < 1) {
    cerr << "Number of rounds must be greater than 0\n";
    return;
  }

  dvectlist vectors;
  long n = dvect_load(in, vectors);

  if (n <= 0) {
    cerr << "Failed to load any input data\n";
    return;
  }

  kmeansresult result(k);
  stats_vector_kmeans(vectors, k, result, rounds);

  cout << "Total data points: " << vectors.size() << endl << endl;

  for (int i = 0; i < k; i++) {
    cout << "Cluster-" << i << ": " << result.get_counts()->at(i) << " entries [ Centroid = "
	 << dvect_tostring(result.get_centroids()->at(i)) << " ]\n";
  }
  cout << endl << "Distortion: " << result.get_distortion() << endl;
  cout << "Average Distortion: " << result.get_distortion() / (long) vectors.size() << endl;
  cout << "Squared Distance Distortion: " << result.get_squared_distance_distortion() << endl;

  if (verbose) {
    cout << endl;
    for (int i = 0; i < k; i++) {
      cout << "Cluster-" << i << ":\n";
      for (long j = 0; j < (long) vectors.size(); j++) {
	if (i == result.get_assignments()->at(j)) {
	  cout << "  " << dvect_tostring(vectors.at(j)) << endl;
	}
      } 
      cout << endl;
    }
  }
}

void compute_vector_xmeans(istream & in, const int k_min, const int k_max, const int rounds, const bool verbose) {
  if (k_min > k_max) {
    cerr << "k_max must be greater than or equal to k_min\n";
    return;
  } else if (k_min < 1) {
    cerr << "k_min must not be greater than 0\n";
  } else if (rounds < 1) {
    cerr << "Number of rounds must be greater than 0\n";
    return;
  }

  dvectlist vectors;
  long n = dvect_load(in, vectors);

  if (n <= 0) {
    cerr << "Failed to load any input data\n";
    return;
  }

  int turns = k_max - k_min + 1;

  kmeansresult* results[turns];
  double bic_scores[turns];

  for (int i = 0; i < turns; i++) {
    int k = k_min + i;
    kmeansresult* current_result = new kmeansresult(k);
    stats_vector_kmeans(vectors, k, *current_result, rounds);
    double bic = stats_vector_kmeans_bic(vectors, *current_result);
    results[i] = current_result;
    bic_scores[i] = bic;
    cout << "k = " << k << "; BIC Score = " << bic << endl;
  }

  double max_bic = -std::numeric_limits<double>::max();
  int max_index = 0;
  for (int i = 0; i < turns; i++) {
    if (bic_scores[i] > max_bic) {
      max_index = i;
      max_bic = bic_scores[i];
    }
  }

  int k = k_min + max_index;

  cout << endl << "Chosen k value: " << k << " (with BIC score: " << max_bic << ")" << endl;
  cout << "Total data points: " << vectors.size() << endl << endl;

  for (int i = 0; i < k; i++) {
    cout << "Cluster-" << i << ": " << (*results[max_index]).get_counts()->at(i) << " entries [ Centroid = "
	 << dvect_tostring((*results[max_index]).get_centroids()->at(i)) << " ]\n";
  }
  cout << endl << "Distortion: " << (*results[max_index]).get_distortion() << endl;
  cout << "Average Distortion: " << (*results[max_index]).get_distortion() / (long) vectors.size() << endl;
  cout << "Squared Distance Distortion: " << (*results[max_index]).get_squared_distance_distortion() << endl;

  if (verbose) {
    cout << endl;
    for (int i = 0; i < k; i++) {
      cout << "Cluster-" << i << ":\n";
      for (long j = 0; j < (long) vectors.size(); j++) {
	if (i == (*results[max_index]).get_assignments()->at(j)) {
	  cout << "  " << dvect_tostring(vectors.at(j)) << endl;
	}
      } 
      cout << endl;
    }
  }

  for (int i = 0; i < turns; i++) {
    delete results[i];
  }
}
