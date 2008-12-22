#ifndef __Polynominal_h_wanghan__
#define __Polynominal_h_wanghan__


#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>

#include <cassert>
#include <cmath>

class Term
{
  std::string name;
  double scale;
  std::vector<int > index;
public:
  Term (const std::string & name, const unsigned & dim);
  Term (const std::string & name, const unsigned & dim, const double & scale);
  Term (const std::string & name, const double & scale, const std::vector<int > & index);
  virtual ~Term () {};
public:
  void reinit (const std::string & name, const double & scale, const std::vector<int > & index);
  bool isZero ();
  void regZero ();

  std::string getName () {return name;};
  double & getScale () {return scale;};
  std::vector<int > & getIndex () {return index;};
  int & getIndex (const unsigned & i) {return index[i];};
  const std::string & getName () const {return name;};
  const double & getScale () const {return scale;};
  const std::vector<int > & getIndex () const {return index;};
  const int & getIndex (const unsigned & i) const {return index[i];};

  void derivative (const unsigned & i);
  Term & operator += (const Term & term);
  Term & operator *= (const Term & term);
  Term & operator *= (const double & scale);

  virtual void print ();
}
    ;

class Polynominal : public std::list<Term >
{
  std::string name;
public:
  Polynominal (const std::string & name);
public:
  Polynominal & operator += (const Term & term);
  Polynominal & operator *= (const Term & term);
  void print ();
  void derivative (const unsigned & i);

  void delZero ();
}
    ;





#endif
