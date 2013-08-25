#ifndef STATS_H_
#define STATS_H_

#include <string>

#include "dvect.h"

double stats_sum(const dvect & v);
double stats_mean(const dvect & v);
double stats_stddev(const dvect & v);
double stats_median(const dvect & v);
double stats_min(const dvect & v);
double stats_max(const dvect & v);

class kmeansresult {
 private:
  lvect counts;
  ivect assignments;
  double sd_distortion;
  double d_distortion;

 public:
  lvect get_counts() { return counts; };
  ivect get_assignments() { return assignments; };
  double get_squared_distance_distortion() { return sd_distortion; }
  double get_distortion() { return d_distortion; }
  void set_counts(lvect c) { counts = c; };
  void set_assignments(ivect a) { assignments = a; };
  void set_squared_distance_distortion(double d) { sd_distortion = d; };
  void set_distortion(double d) { d_distortion = d; };
};

class kmeansresult2 {
 private:
  lvect *counts;
  ivect *assignments;
  dvectlist *centroids;
  double sd_distortion;
  double d_distortion;
  bool assignments_set;
  bool centroids_set;

 public:
  kmeansresult2(int k);
  ~kmeansresult2();
  const lvect* const get_counts() { return counts; };
  const ivect* const get_assignments() { return assignments; };
  const dvectlist* const get_centroids() { return centroids; };
  const double get_squared_distance_distortion() { return sd_distortion; };
  const double get_distortion() { return d_distortion; };
  void set_assignments(ivect* a);
  void set_centroids(dvectlist* c);
  void set_squared_distance_distortion(double d) { sd_distortion = d; };
  void set_distortion(double d) { d_distortion = d; };
  void increment_count(int index);
};

dvect stats_vector_centroid(const dvectlist & vectors);
dvectlist stats_vector_kmeans(const dvectlist & vectors, const int k, kmeansresult & result, const int rounds = 20);

#endif
