#ifndef __LQ_Factoriztion_wanghan__
#define __LQ_Factoriztion_wanghan__

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "mkl.h"

template <typename ValueType>
class LQFactorization
{
public:
  //void reinit ();
  int run (const std::vector<ValueType > & matrix, 
	   const int & m, 
	   const int & n,
	   std::vector<ValueType > & L,
	   std::vector<ValueType > & Q);
}
    ;


////////////////////////////////////////
// specify template
template <typename ValueType>
int LQFactorization<ValueType >::run (const std::vector<ValueType > & matrix, 
				      const int & m, 
				      const int & n,
				      std::vector<ValueType > & L,
				      std::vector<ValueType > & Q)
{
  std::cout << "LQFactorization: the template parameter has not been implemented, do nothing\n" ;
  return -1;
}

template <>
class LQFactorization<double > 
{
public:
  //void reinit ();
  int run (const std::vector<double > & matrix, 
	   const int & m, 
	   const int & n,
	   std::vector<double > & L,
	   std::vector<double > & Q);
}
    ;









#endif
