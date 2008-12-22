#ifndef __LU_solver_h_wanghan__
#define __LU_solver_h_wanghan__

#include"mkl.h"

class LUsolver
{
  std::vector<double > a;
  std::vector<int >    ipiv;
  int dim;
  int nrhs;
  char trans;
public:
  void reinit (const std::vector<double > & matrix, const int & dim);
  void solve  (std::vector<double > & result);
  void solve  (const std::vector<double > & rhs, std::vector<double > & result);
}
    ;

#endif


void LUsolver::reinit (const std::vector<double > & matrix, const int & dim_)
{
  dim = dim_;
  nrhs = 1;
  trans = 'T';
  
  a.resize (dim*dim);
  a = matrix;
  ipiv.resize (dim);

  int info;
  
  if (matrix.size() != unsigned (dim * dim)){
    std::cerr << "the dimension of the matrix is wrong!\n";
    exit(1);
  }
  if (dim_ < 0){
    std::cerr << "dimension should not be negative!\n";
    exit (1);
  }

  dgetrf (&dim, &dim, &a[0], &dim, &ipiv[0], &info);
}

void LUsolver::solve (std::vector<double > & result)
{
  result.resize (dim);

  int info;
  
  dgetrs (&trans, &dim, &nrhs, &a[0], &dim, &ipiv[0], &result[0], &dim, &info);
}

void LUsolver::solve (const std::vector<double > & rhs, std::vector<double > & result)
{
  result = rhs;
  if (result.size() != unsigned(dim)){
    std::cerr << "the dimension of rhs is wrong!\n";
    exit(1);
  }
  
  int info;
  
  dgetrs (&trans, &dim, &nrhs, &a[0], &dim, &ipiv[0], &result[0], &dim, &info);
}

  
  

