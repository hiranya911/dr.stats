#ifndef DVECT_H_
#define DVECT_H_

#include <iostream>
#include <vector>

typedef std::vector< double > dvect;
typedef std::vector< dvect > dvectlist;
typedef std::vector< double >::iterator dvectiter;
typedef std::vector< double >::const_iterator dvectciter;
typedef std::vector< dvect >::iterator dvectlistiter;
typedef std::vector< dvect >::const_iterator dvectlistciter;
typedef std::vector< long > lvect;
typedef std::vector< long >::iterator lvectiter;
typedef std::vector< int > ivect;
typedef std::vector< int >::iterator ivectiter;

void dvect_print(const dvect & v);
std::string dvect_tostring(const dvect & v);
long dvect_load(std::istream & in, dvect & numbers);
long dvect_load(std::istream & in, dvectlist & vectors);
long dvect_load(const std::string & str, dvect & numbers);

#endif
