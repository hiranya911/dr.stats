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
double within_cluster_squared_sum(const dvectlist & vectors, const ivect* const assignments, const dvectlist* const centroids);
double within_cluster_distance(const dvectlist & vectors, const ivect* const assignments, const dvectlist* const centroids);
void update_centroids(dvectlist* const centroids, const dvectlist & vectors, const ivect* const assignments);

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
  double max = -std::numeric_limits<double>::max();
  for (dvectciter it = v.begin(); it != v.end(); it++) {
    if (*it > max) {
      max = *it;
    }
  }
  return max;
}

void stats_sort(dvect & v) {
  std::sort(v.begin(), v.end());
}

dvect stats_vector_add(const dvectlist & vectors) {
  long dvect_size = (long) vectors.front().size();
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
  }
  return result;
}

dvect stats_vector_centroid(const dvectlist & vectors) {
  long size = (long) vectors.size();
  long dvect_size = (long) vectors.front().size();
  dvect result = stats_vector_add(vectors);
  for (long i = 0; i < dvect_size; i++) {
    result[i] /= size;
  }
  return result;
}

double stats_vector_dot_product(const dvect & v1, const dvect & v2) {
  if (v1.size() != v2.size()) {
    throw std::runtime_error("Vector sizes do not match");
  }
  long size = (long) v1.size();
  double sum = 0.0;
  for (long i = 0; i < size; i++) {
    sum += v1.at(i) * v2.at(i);
  }
  return sum;
} 

void stats_vector_kmeans(const dvectlist & vectors, const int k, kmeansresult & result, const int rounds) {
  if (k < 1) {
    throw std::runtime_error("k must not be less than 1");
  }

  long dvect_size = (long) vectors.front().size();
  long n = (long) vectors.size();

  if (k > n) {
    throw std::runtime_error("k must not be greater than the number of data points");
  }
  
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
      ivect* assignments = new ivect(n, -1);
      // Assign each vector to a cluster
      for (long v = 0; v < n; v++) {
	long cluster = 0;
	double min_squared_distance = std::numeric_limits<double>::max();
	for (int i = 0; i < k; i++) {
	  double distance = stats_vector_euclidean_distance(vectors.at(v), centroids->at(i));
	  double squared_distance = distance * distance;
	  if (squared_distance < min_squared_distance) {
	    cluster = i;
	    min_squared_distance = squared_distance;
	  }
	}
	(*assignments)[v] = cluster;
      }
      
      double new_wcss = within_cluster_squared_sum(vectors, assignments, centroids);
      if (new_wcss != wcss) {
	// WCSS has changed since last iteration - Need to run another iteration
	wcss = new_wcss;
	update_centroids(centroids, vectors, assignments);
	delete assignments;
      } else {
	// Algorithm has converged to a solution
	if (wcss < result.get_squared_distance_distortion()) {
	  // We have a better solution than the previous rounds - Store the latest result
	  result.set_squared_distance_distortion(wcss);
	  result.set_centroids(centroids);
	  result.set_assignments(assignments);
	} else {
	  // We already have a better solution from a previous round - Ignore the latest result
	  delete assignments;
	  delete centroids;
	}
	break;
      }
    }
  }

  double distortion = within_cluster_distance(vectors, result.get_assignments(), result.get_centroids());
  for (long i = 0; i < n; i++) {
    result.increment_count(result.get_assignments()->at(i));
  }
  result.set_distortion(distortion);
}

double stats_vector_kmeans_bic(const dvectlist & vectors, const kmeansresult & result) {
  long n = (long) vectors.size();
  long m = (long) vectors.front().size();
  int k = (int) result.get_centroids()->size();

  if (n != (long) result.get_assignments()->size()) {
    std::stringstream ss;
    ss << "Assignment count mismatch; Expected: " << n << ", Found: " << result.get_assignments()->size();
    throw std::runtime_error(ss.str());
  }

  const double ROOT_2_PI = sqrt(2 * 3.141592653589793);

  double sigma_squared = (result.get_squared_distance_distortion()/(n - k));
  double sigma = sqrt(sigma_squared);

  double log_likelihood = 0.0;
  double A = log(1/(ROOT_2_PI * pow(sigma, m)));

  for (int i = 0; i < n; i++) {
    dvect vector = vectors.at(i);
    if (m != (long) vector.size()) {
      std::stringstream ss;
      ss << "Vector size mismatch; Expected: " << m << ", Found: " << vector.size();
      throw std::runtime_error(ss.str());
    }
    int cluster = result.get_assignments()->at(i);
    dvect centroid = result.get_centroids()->at(cluster);
    double distance = stats_vector_euclidean_distance(vector, centroid);
    double B = (1/(2 * sigma_squared)) *  distance * distance;
    double C = log((result.get_counts()->at(cluster))/(double) n);
    log_likelihood += (A - B + C);
  }

  long p = k + m * k;
  double BIC = log_likelihood - ((p/2) * log(n));
  return BIC;
}

double stats_vector_euclidean_distance(const dvect & v1, const dvect & v2) {
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

void update_centroids(dvectlist* const centroids, const dvectlist & vectors, const ivect* const assignments) {
  int k = (int) centroids->size();
  long v_size = (long) vectors.front().size();
  long *counts = new long[k];
  for (int i = 0; i < k; i++) {
    counts[i] = 0L;
    (*centroids)[i] = dvect(v_size, 0.0);
  }
  for (long i = 0; i < (long) vectors.size(); i++) {
    int cluster = assignments->at(i);
    for (long j = 0; j < v_size; j++) {
      centroids->at(cluster)[j] += vectors.at(i).at(j);
    }
    counts[cluster]++;
  }
  for (int i = 0; i < k; i++) {
    for (long j = 0; j < v_size; j++) {
      centroids->at(i)[j] /= counts[i];
    }
  }
  delete [] counts;
}

double within_cluster_squared_sum(const dvectlist & vectors, const ivect* const assignments, const dvectlist* const centroids) {
  double wcss = 0.0;
  for (long i = 0; i < (long) vectors.size(); i++) {
    double distance = stats_vector_euclidean_distance(vectors.at(i), centroids->at(assignments->at(i)));
    wcss += distance * distance;
  }
  return wcss;
}

double within_cluster_distance(const dvectlist & vectors, const ivect* const assignments, const dvectlist* const centroids) {
  double wcd = 0.0;
  for (long i = 0; i < (long) vectors.size(); i++) {
    wcd += stats_vector_euclidean_distance(vectors.at(i), centroids->at(assignments->at(i)));
  }
  return wcd;
}

void validate_vector(const dvect & v) {
  if (v.empty()) {
    throw std::runtime_error("Empty input vector");
  }
}
