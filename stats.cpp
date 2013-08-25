#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <stdexcept>

#include "dvect.h"
#include "stats.h"

using std::string;

void validate_vector(const dvect & v);
double within_cluster_squared_sum(const dvectlist & vectors, const ivect & assignments, const dvectlist & centroids);
double within_cluster_distance(const dvectlist & vectors, const ivect & assignments, const dvectlist & centroids);
void update_centroids(dvectlist & centroids, const dvectlist & vectors, const ivect & assignments);
double euclidean_distance(const dvect & v1, const dvect & v2);

kmeansresult::kmeansresult(int k) {
  counts = new lvect(k, 0L);
  assignments_set = false;
  centroids_set = false;
}

kmeansresult::~kmeansresult() {
  delete counts;
  if (assignments_set) {
    delete assignments;
  }
  if (centroids_set) {
    delete centroids;
  }
}

void kmeansresult::set_assignments(ivect* a) {
  if (assignments_set) {
    delete assignments;
  }
  assignments = a;
  assignments_set = true;
}

void kmeansresult::set_centroids(dvectlist* c) {
  if (centroids_set) {
    delete centroids;
  }
  centroids = c;
  centroids_set = true;
}

void kmeansresult::increment_count(int index) {
  counts->at(index)++;
}

double stats_sum(const dvect & v) {
  double sum = 0.0;
  for (dvectciter it = v.begin(); it != v.end(); it++) {
    sum += (*it);
  }
  return sum;
}

double stats_mean(const dvect & v) {
  validate_vector(v);
  return stats_sum(v) / (long) v.size();
}

double stats_stddev(const dvect & v) {
  validate_vector(v);
  double ss = 0.0;
  double mean = stats_mean(v);
  for (dvectciter it = v.begin(); it != v.end(); it++) {
    ss += (*it - mean) * (*it - mean);
  }
  return sqrt(ss / (long) v.size());
}

double stats_median(const dvect & v) {
  validate_vector(v);
  dvect numbers = v;
  std::sort(numbers.begin(), numbers.end());
  double median;
  long size = (long) v.size();
  if (size % 2 == 0) {
    double med1 = numbers.at(size / 2);
    double med2 = numbers.at((size / 2) - 1);
    return (med1 + med2) / 2.0;
  } else {
    return numbers.at(floor(size / 2.0));
  }
}

double stats_min(const dvect & v) {
  validate_vector(v);
  double min = std::numeric_limits<double>::max();
  for (dvectciter it = v.begin(); it != v.end(); it++) {
    if (*it < min) {
      min = *it;
    }
  }
  return min;
}

double stats_max(const dvect & v) {
  validate_vector(v);
  double max = std::numeric_limits<double>::min();
  for (dvectciter it = v.begin(); it != v.end(); it++) {
    if (*it > max) {
      max = *it;
    }
  }
  return max;
}

dvect stats_vector_centroid(const dvectlist & vectors) {
  long size = (long) vectors.size();
  long dvect_size = (long) vectors.front().size();
  long vcount = 1L;
  dvect result(dvect_size, 0.0);
  for (dvectlistciter it = vectors.begin(); it != vectors.end(); it++) {
    if (dvect_size != (long) (*it).size()) {
      std::stringstream ss;
      ss << "Vector size mismatch; Expected: " << dvect_size << ", Found: " << (*it).size();
      throw std::runtime_error(ss.str());
    }
    for (long i = 0; i < dvect_size; i++) {
      result[i] += (*it)[i];
    }
    vcount++;
  }
  for (long i = 0; i < dvect_size; i++) {
    result[i] /= size;
  }
  return result;
}

void stats_vector_kmeans(const dvectlist & vectors, const int k, kmeansresult & result, const int rounds) {
  long dvect_size = (long) vectors.front().size();
  long n = (long) vectors.size();
  
  result.set_squared_distance_distortion(std::numeric_limits<double>::max());
  srand(time(NULL));

  for (int rounds_count = 0; rounds_count < rounds; rounds_count++) {
    dvectlist* centroids = new dvectlist;
    centroids->reserve(k);

    // Initialize the clusters
    std::set<int> centroid_indices;
    while ((int) centroid_indices.size() < k) {
      centroid_indices.insert(rand() % n);
    }
    for (std::set<int>::iterator it = centroid_indices.begin(); it != centroid_indices.end(); it++) {
      centroids->push_back(vectors.at(*it));
    }

    double wcss = -1.0;
    while (true) {
      ivect* assignments = new ivect;
      assignments->reserve(n);
      for (dvectlistciter it = vectors.begin(); it != vectors.end(); it++) {
	long cluster = 0;
	double min_squared_distance = std::numeric_limits<double>::max();
	for (int i = 0; i < k; i++) {
	  double distance = euclidean_distance(*it, centroids->at(i));
	  double squared_distance = distance * distance;
	  if (squared_distance < min_squared_distance) {
	    cluster = i;
	    min_squared_distance = squared_distance;
	  }
	}
	assignments->push_back(cluster);
      }
      
      double new_wcss = within_cluster_squared_sum(vectors, *assignments, *centroids);
      if (new_wcss != wcss) {
	wcss = new_wcss;
	update_centroids(*centroids, vectors, *assignments);
	delete assignments;
      } else {
	if (wcss < result.get_squared_distance_distortion()) {
	  result.set_squared_distance_distortion(wcss);
	  result.set_centroids(centroids);
	  result.set_assignments(assignments);
	} else {
	  delete assignments;
	  delete centroids;
	}
	break;
      }
    }
  }

  double distortion = within_cluster_distance(vectors, *(result.get_assignments()), *(result.get_centroids()));
  for (long i = 0; i < n; i++) {
    result.increment_count(result.get_assignments()->at(i));
  }
  result.set_distortion(distortion);
}

void update_centroids(dvectlist & centroids, const dvectlist & vectors, const ivect & assignments) {
  int k = (int) centroids.size();
  long v_size = (long) vectors.front().size();
  long *counts = new long[k];
  for (int i = 0; i < k; i++) {
    counts[i] = 0L;
    centroids[i] = dvect(v_size, 0.0);
  }
  for (long i = 0; i < (long) vectors.size(); i++) {
    int cluster = assignments.at(i);
    for (long j = 0; j < v_size; j++) {
      centroids.at(cluster)[j] += vectors.at(i).at(j);
    }
    counts[cluster]++;
  }
  for (int i = 0; i < k; i++) {
    for (long j = 0; j < v_size; j++) {
      centroids.at(i)[j] /= counts[i];
    }
  }
  delete [] counts;
}

double within_cluster_squared_sum(const dvectlist & vectors, const ivect & assignments, const dvectlist & centroids) {
  double wcss = 0.0;
  for (long i = 0; i < (long) vectors.size(); i++) {
    double distance = euclidean_distance(vectors.at(i), centroids.at(assignments.at(i)));
    wcss += distance * distance;
  }
  return wcss;
}

double within_cluster_distance(const dvectlist & vectors, const ivect & assignments, const dvectlist & centroids) {
  double wcd = 0.0;
  for (long i = 0; i < (long) vectors.size(); i++) {
    wcd += euclidean_distance(vectors.at(i), centroids.at(assignments.at(i)));
  }
  return wcd;
}

double euclidean_distance(const dvect & v1, const dvect & v2) {
  if (v1.size() != v2.size()) {
    std::stringstream ss;
    ss << "Vector size mismatch; Expected: " << v1.size() << ", Found: " << v2.size();
    throw std::runtime_error(ss.str());
  }
  double ss = 0.0;
  for (long i = 0; i < (long) v1.size(); i++) {
    ss += (v1.at(i) - v2.at(i)) * (v1.at(i) - v2.at(i));
  }
  return sqrt(ss); 
}

void validate_vector(const dvect & v) {
  if (v.empty()) {
    throw std::runtime_error("Empty input vector");
  }
}
