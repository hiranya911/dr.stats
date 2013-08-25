#ifndef DVECT_H_
#define DVECT_H_

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

void dvect_init(dvect & v, const std::string & line);
void dvect_print(const dvect & v);
std::string dvect_tostring(const dvect & v);
long dvect_load(const std::string & file, dvect & numbers);
long dvect_load(const std::string & file, dvectlist & vectors);

#endif
