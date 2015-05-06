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
