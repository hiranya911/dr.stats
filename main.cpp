#include <iostream>

#include "boost/program_options.hpp"
#include "dvect.h"
#include "stats.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

void compute_simple_mode(const string & file);
void compute_vector_centroid(const string & file);
void compute_vector_kmeans(const string & file, const int k);

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
    ("input,i", po::value<string>(&input_file)->required(), "Input file")
    ("clusters,k", po::value<int>(), "Number of clusters (used for kmeans)")
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
    if (mode == "simple") {
      compute_simple_mode(input_file);
    } else if (mode == "cent") {
      compute_vector_centroid(input_file);
    } else if (mode == "kmeans") {
      if (vm.count("clusters")) {
	int k = vm["clusters"].as<int>();
	compute_vector_kmeans(input_file, k);
      } else {
	cerr << "Error: k is required but missing\n";
	return 1;
      } 
    } else {
      cerr << "Unsupported calculation mode: " << mode << endl;
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

void compute_simple_mode(const string & file) {
  dvect numbers;
  long size = dvect_load(file, numbers);

  if (size == 0) {
    cout << "No data available in file: " << file << endl;
    return;
  } else if (size < 0) {
    cout << "Failed to load the file: " << file << endl;
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

void compute_vector_centroid(const string & file) {
  dvectlist vectors;
  long size = dvect_load(file, vectors);

  if (size == 0) {
    cout << "No data available in file: " << file << endl;
    return;
  } else if (size < 0) {
    cout << "Failed to load the file: " << file << endl;
    return;
  }

  dvect result = stats_vector_centroid(vectors);
  cout << "Number of vectors: " << size << endl;
  cout << "Vector dimensions: " << vectors.front().size() << endl;
  cout << "Centroid: ";
  dvect_print(result);
}

void compute_vector_kmeans(const string & file, const int k) {
  if (k < 1) {
    cerr << "Cluster count (k) must be greater than 0\n";
    return;
  }
  dvectlist vectors;
  long n = dvect_load(file, vectors);

  if (n == 0) {
    cout << "No data available in file: " << file << endl;
    return;
  } else if (n < 0) {
    cout << "Failed to load the file: " << file << endl;
    return;
  }

  kmeansresult result(k);
  stats_vector_kmeans(vectors, k, result);

  for (int i = 0; i < k; i++) {
    cout << "Cluster-" << i << ": " << result.get_counts()->at(i) << " entries [ Centroid = "
	 << dvect_tostring(result.get_centroids()->at(i)) << " ]\n";
  }
  cout << endl << "Distortion: " << result.get_distortion() << endl;
  cout << "Squared Distance Distortion: " << result.get_squared_distance_distortion() << endl;
}
