#ifndef __VectorOperation_h_wanghan__
#define __VectorOperation_h_wanghan__

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


namespace VectorOperation{

// level 1
    // set all v component 0
    void zero (std::vector<double > & v);
    // v = alpha * v
    void scale (std::vector<double > & v, const double & alpha);
    // v += alpha * a
    void add (std::vector<double > & v, 
	      const double & alpha, const std::vector<double > & a);
    // v += alpha * a + beta * b
    void add (std::vector<double > & v,
	      const double & alpha, const std::vector<double > & a,
	      const double & beta,  const std::vector<double > & b);
    // v = s * v + alpha * a
    void sadd (std::vector<double > & v, const double & s, 
	       const double & alpha, const std::vector<double > & a);
    // v = s * v + alpha * a + beta * b
    void sadd (std::vector<double > & v, const double & s,
	       const double & alpha, const std::vector<double > & a,
	       const double & beta,  const std::vector<double > & b);
    // return u .* v
    double dot (const std::vector<double > & u, const std::vector<double > & v);
    // print to stdout
    void print (const std::vector<double > & v, 
		const bool & inColume = true);
    // print to file
    void print (const std::string & filename, const std::vector<double > & v);
    // read from a file
    bool read (const std::string & filename, std::vector<double > & v);

// levle 2
    // v = beta * v + alpha * A^(t) * a
    void mv (std::vector<double > & v,
	     const double & beta,
	     const double & alpha, 
	     const std::vector<double > & A, const bool & transA, 
	     const std::vector<double > & a);
    // v = alpha * A^(t) * a
    void mv (std::vector<double > & v,
	     const double & alpha, 
	     const std::vector<double > & A, const bool & transA, 
	     const std::vector<double > & a);
    // return u^(t) * A * v
    double quadForm (const std::vector<double > & u, const std::vector<double > & A, std::vector<double > & v);
// level 3
    // cal the dimension of A
    unsigned dim (const std::vector<double > & A, bool & dim_valid);
    // A += alpha * I
    void eyeCorrect (std::vector<double > & A, const double & alpha);
    // for multiply of two square matrices, C = alpha * A^(t) * B^(t)
    void mm (std::vector<double > & C, 
	     const double & alpha, 
	     const std::vector<double > & A, const bool & tranA, 
	     const std::vector<double > & B, const bool & tranB);
// printers and loaders
    // print matrix 
    void printMatrix (const std::vector<double > & A);
    void printMatrix (const std::vector<double > & A,
		      const int & m, 
		      const int & n );
    void printMatrix (const std::string & filename, 
		      const std::vector<double > & A);
    void printMatrix (const std::string & filename, 
		      const std::vector<double > & A,
		      const int & m, 
		      const int & n );
    void printVector (const std::vector<double > & v);
    void printVector (const std::vector<double > & u, 
		      const std::vector<double > & v);
    void printVector (const std::string & filename, 
		      const std::vector<double > & v);
    void printVector (const std::string & filename, 
		      const std::vector<double > & u, 
		      const std::vector<double > & v);
    // load vector
    void loadVector (const std::string & filename, 
		     std::vector<double > & v );
    void loadVector (const std::string & filename, 
		     std::vector<double > & u ,
		     std::vector<double > & v );
};



#endif

