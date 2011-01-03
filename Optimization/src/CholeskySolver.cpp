#include "CholeskySolver.h"

int CholeskySolver<double >::reinit (const std::vector<double > & matrix, const int & dim_)
{
  uplo = 'U';
  dim = dim_;
  nrhs = 1;
  
  if (matrix.size() != unsigned (dim * dim)){
    std::cerr << "CholeskySolver: the dimension of the matrix is wrong!\n";
    exit(1);
  }
  if (dim_ < 0){
    std::cerr << "CholeskySolver: dimension should not be negative!\n";
    exit (1);
  }

  a = matrix;
  
  int info;
  
  dpotrf (&uplo, &dim, &a[0], &dim, &info);

  return info;
}

int CholeskySolver<double >::solve (std::vector<double > &result)
{
  result.resize (dim);

  int info;
  
  dpotrs (&uplo, &dim, &nrhs, &a[0], &dim, &result[0], &dim, &info);

  return info;
}

int CholeskySolver<double >::solve (const std::vector<double > & rhs, std::vector<double > & result)
{  
  result = rhs;
  
  if (result.size() != unsigned(dim)){
    std::cerr << "CholeskySolver: the dimension of rhs is wrong!\n";
    exit(1);
  }

  int info;

  dpotrs (&uplo, &dim, &nrhs, &a[0], &dim, &result[0], &dim, &info);

  return info;
}

