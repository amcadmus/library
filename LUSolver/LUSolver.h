#ifndef __LU_Solver_h_wanghan__
#define __LU_Solver_h_wanghan__

#include "mkl.h"
#include <iostream>
#include <vector>

template <typename ValueType>
class LUSolver
{
  std::vector<ValueType > a;
  std::vector<int >    ipiv;
  int dim;
  int nrhs;
  char trans;
public:
  int reinit (const std::vector<ValueType > & matrix, const int & dim);
  int solve  (std::vector<ValueType > & result);
  int solve  (const std::vector<ValueType > & rhs, std::vector<ValueType > & result);
}
    ;

template <typename ValueType>
int LUSolver<ValueType >::reinit (const std::vector<ValueType > & matrix, const int & dim_)
{
  std::cout << "LUSolver: the template parameter is not realized, do nothing!\n";
  return -1;
}

template <typename ValueType>
int LUSolver<ValueType >::solve (std::vector<ValueType > & result)
{
  std::cout << "LUSolver: the template parameter is not realized, do nothing!\n";
  return -1;
}

template <typename ValueType>
int LUSolver<ValueType >::solve (const std::vector<ValueType > & rhs, std::vector<ValueType > & result)
{
  std::cout << "LUSolver: the template parameter is not realized, do nothing!\n";
  return -1;
}


  
template <>
class LUSolver<double >
{
  std::vector<double > a;
  std::vector<int >    ipiv;
  int dim;
  int nrhs;
  char trans;
public:
  int reinit (const std::vector<double > & matrix, const int & dim);
  int solve  (std::vector<double > & result);
  int solve  (const std::vector<double > & rhs, std::vector<double > & result);
}
    ;

  
#endif


