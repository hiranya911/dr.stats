#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "dvect.h"

using namespace std;

dvect centroid(const dvectlist & v);
double euclidean_distance(const dvect & v1, const dvect & v2);

int main() {
  cout.precision(10);
  ifstream file("example.txt", ios::in);
  if (file) {
    string line;
    string buffer;
    dvectlist list;
    while (getline(file, line)) {
      dvect elements;
      dvect_init(elements, line);
      list.push_back(elements);
      dvect_print(elements);
    }
    file.close();
    cout << "Average: ";
    dvect_print(centroid(list));

    dvect v1, v2;
    dvect_init(v1, "3 4");
    dvect_init(v2, "0 0");
    cout << "Distance: " << euclidean_distance(v1, v2) << endl;
  } else {
    cerr << "Failed to open file!\n";
  }
  return 0;
}

dvect centroid(const dvectlist & v) {
  int size = v.at(0).size();
  for (int i = 1; i < v.size(); i++) {
    if (size != v.at(i).size()) {
      throw "Vector size mismatch";
    }
  }

  dvect result(size, 0.0);
  for (int i = 0; i < v.size(); i++) {
    for (int j = 0; j < size; j++) {
      result[j] += v.at(i).at(j);
    }
  }
  for (int i = 0; i < size; i++) {
    result[i] /= v.size();
  }
  return result;
}

double euclidean_distance(const dvect & v1, const dvect & v2) {
  if (v1.size() != v2.size()) {
    throw "Vector size mismatch";
  }
  double ss = 0.0;
  for (int i = 0; i < v1.size(); i++) {
    ss += (v1.at(i) - v2.at(i)) * (v1.at(i) - v2.at(i));
  }
  return sqrt(ss);
}

