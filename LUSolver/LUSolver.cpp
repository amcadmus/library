#include "LUSolver.h"

int LUSolver<double >::reinit (const std::vector<double > & matrix, const int & dim_)
{
  dim = dim_;
  nrhs = 1;
  trans = 'T';

  if (matrix.size() != unsigned (dim * dim)){
    std::cerr << "LUSolver: the dimension of the matrix is wrong!\n";
    exit(1);
  }
  if (dim_ < 0){
    std::cerr << "LUSolver: dimension should not be negative!\n";
    exit (1);
  }
  
  a.resize (dim*dim);
  a = matrix;
  ipiv.resize (dim);

  int info;

  dgetrf (&dim, &dim, &a[0], &dim, &ipiv[0], &info);

  return info;
}

int LUSolver<double >::solve (std::vector<double > & result)
{
  result.resize (dim);

  int info;
  
  dgetrs (&trans, &dim, &nrhs, &a[0], &dim, &ipiv[0], &result[0], &dim, &info);

  return info;
}

int LUSolver<double >::solve (const std::vector<double > & rhs, std::vector<double > & result)
{
  result = rhs;
  if (result.size() != unsigned(dim)){
    std::cerr << "LUSolver: the dimension of rhs is wrong!\n";
    exit(1);
  }
  
  int info;
  
  dgetrs (&trans, &dim, &nrhs, &a[0], &dim, &ipiv[0], &result[0], &dim, &info);

  return info;
}

  
