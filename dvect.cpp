#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "dvect.h"

using std::cout;
using std::string;

void dvect_init(dvect & v, const string & line) {
  std::stringstream ss(line);
  double d;
  while (ss >> d) {
    v.push_back(d);
  }
}

void dvect_print(const dvect & v) {
  cout << dvect_tostring(v) << std::endl;
}

string dvect_tostring(const dvect & v) {
  std::stringstream ss;
  ss << "(";
  for (dvectciter it = v.begin(); it != v.end(); ++it) {
    if (it != v.begin()) {
      ss << ",";
    }
    ss << *it;
  }
  ss << ")";
  return ss.str();
}

long dvect_load(const string & file, dvect & numbers) {
  long count = 0L;
  std::ifstream input_file(file.c_str(), std::ifstream::in);
  if (input_file) {
    double n;
    while (input_file >> n) {
      numbers.push_back(n);
      count++;
    }
    input_file.close();
  } else {    
    count = -1;
  }
  return count;
}

long dvect_load(const string & file, dvectlist & vectors) {
  long count = 0L;
  std::ifstream input_file(file.c_str(), std::ifstream::in);
  if (input_file) {
    string line;
    while (getline(input_file, line)) {
      if (line.empty()) {
	continue;
      }
      dvect v;
      dvect_init(v, line);
      vectors.push_back(v);
      count++;
    }
    input_file.close();
  } else {
    count = -1;
  }
  return count;
}
