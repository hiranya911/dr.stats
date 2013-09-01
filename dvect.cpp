#include <iostream>
#include <sstream>
#include <vector>

#include "dvect.h"

void dvect_print(const dvect & v) {
  std::cout << dvect_tostring(v) << std::endl;
}

std::string dvect_tostring(const dvect & v) {
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

long dvect_load(std::istream & in, dvect & numbers) {
  long count = 0L;
  double n;
  while (in >> n) {
    numbers.push_back(n);
    count++;
  }
  return count;
}

long dvect_load(const std::string & str, dvect & numbers) {
  std::stringstream ss(str);
  return dvect_load(ss, numbers);
}

long dvect_load(std::istream & in, dvectlist & vectors) {
  long count = 0L;
  std::string line;
  while (getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    dvect v;
    dvect_load(line, v);
    vectors.push_back(v);
    count++;
  }
  return count;
}
