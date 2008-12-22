#ifndef __Matrix_h_wanghan__
#define __Matrix_h_wanghan__

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
#include <numeric>

#include <cmath>
#include <cassert> 

#include "mkl.h"


class Matrix : public std::vector<double > 
{
private:
  unsigned dim0;
  unsigned dim1;
public:
  Matrix ();
  Matrix (const unsigned & dim1, const unsigned & dim0);
  Matrix (const unsigned & dim1, const unsigned & dim0, const double & value);
  Matrix (const Matrix & mat);
public:
  void reinit (const unsigned & dim1, const unsigned & dim0);
  void reinit (const unsigned & dim1, const unsigned & dim0, const double & value);
public:
  const unsigned & d0 () const {return dim0;}
  const unsigned & d1 () const {return dim1;}
  unsigned & d0 () {return dim0;}
  unsigned & d1 () {return dim1;}
public:
  const double & operator () (const unsigned & i, const unsigned & j) const {
    return this->operator[] (i * dim0 + j);
  }
  double & operator () (const unsigned & i, const unsigned & j) {
    return this->operator[] (i * dim0 + j);
  }
}
    ;

namespace MatrixOperation {
// level 2
    // v = beta * v + alpha * A^(t) * a
    void mv (std::vector<double > & v,
	     const double & beta,
	     const double & alpha,
	     const Matrix & A, 
	     const bool & trans,
	     const std::vector<double > a);
    // v = alpha * A^(t) * a
    void mv (std::vector<double > & v,
	     const double & alpha,
	     const Matrix & A, 
	     const bool & trans,
	     const std::vector<double > a);
    // return u^(t) * A * v
    double quadForm (const std::vector<double > & u,
		     const Matrix & A,
		     const std::vector<double > & v);
// level 3
    // C = alpha * A^(t) * B^(t)
    void mm (Matrix & C,
	     const double & alpha,
	     const Matrix & A, const bool & transA,
	     const Matrix & B, const bool & transB);
    
    void print (const Matrix & mat);
    void print (const std::string & filename, const Matrix & A);
};   
  

#endif
