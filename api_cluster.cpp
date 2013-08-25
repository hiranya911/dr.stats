#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>

#include "stats.h"

using namespace std;

int main(int argc, char** argv) {
  ifstream file(argv[1], ifstream::in);
  int k = atoi(argv[2]);
  if (file) {
    dvectlist vectors;
    vector<string> names;
    double value;
    string name;
    while (file >> value >> name) {
      dvect vector(1, value);
      vectors.push_back(vector);
      names.push_back(name);
    }
    file.close();

    kmeansresult result(k);
    stats_vector_kmeans(vectors, k, result, 1000);
    for (int i = 0; i < k; i++) {
      cout << "Cluster-" << i << ": " << result.get_counts()->at(i) << " entries [ Centroid = "
	   << dvect_tostring(result.get_centroids()->at(i)) << " ]\n";
    }

    cout << "Average Distortion: " << result.get_distortion() / (long) vectors.size() << endl;    
    cout << endl;
    for (int i = 0; i < k; i++) {
      cout << "Cluster-" << i << ":\n";
      for (long j = 0; j < (long) vectors.size(); j++) {
	if (i == result.get_assignments()->at(j)) {
	  cout << "  " << dvect_tostring(vectors.at(j)) << " " << names.at(j) << endl;
	}
      } 
      cout << endl;
    }

  } else {
    cerr << "Failed to open the file: " << argv[1] << endl;
    return 1;
  }
  return 0;
}
