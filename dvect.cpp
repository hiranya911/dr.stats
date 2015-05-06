/*
 * Copyright 2015 Hiranya Jayathilaka
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 *limitations under the License.
 */

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
