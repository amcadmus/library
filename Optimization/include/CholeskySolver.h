#ifndef __Cholesky_solver_h_wanghan__
#define __Cholesky_solver_h_wanghan__

#include"mkl.h"
#include<vector>
#include<iostream>

template <typename ValueType>
class CholeskySolver
{
  std::vector<ValueType > a;
  int dim;
  int nrhs;
  char uplo;
public:
  int reinit (const std::vector<ValueType > & matrix, const int & dim);
  int solve  (std::vector<ValueType > & result);
  int solve  (const std::vector<ValueType > & rhs, std::vector<ValueType > & result);
}
    ;

template <typename ValueType>
int CholeskySolver<ValueType >::reinit (const std::vector<ValueType > & matrix, const int & dim_)
{
  std::cout << "CholeskySolver: the template parameter is not realized, do nothing!\n";
  return -1;
}

template <typename ValueType>
int CholeskySolver<ValueType >::solve (std::vector<ValueType > &result)
{
  std::cout << "CholeskySolver: the template parameter is not realized, do nothing!\n";
  return -1;
}

template <typename ValueType>
int CholeskySolver<ValueType >::solve (const std::vector<ValueType > & rhs, std::vector<ValueType > & result)
{  
  std::cout << "CholeskySolver: the template parameter is not realized, do nothing!\n";
  return -1;
}


template <>
class CholeskySolver<double >
{
  std::vector<double > a;
  int dim;
  int nrhs;
  char uplo;
public:
  int reinit (const std::vector<double > & matrix, const int & dim);
  int solve  (std::vector<double > & result);
  int solve  (const std::vector<double > & rhs, std::vector<double > & result);
}
    ;



#endif
