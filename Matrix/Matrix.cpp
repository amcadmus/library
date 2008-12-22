#include "Matrix.h"

Matrix::Matrix ()
    : dim0 (0), dim1(0)
{
}

Matrix::Matrix (const unsigned & dim1_, const unsigned & dim0_)
    : std::vector<double > (dim0_ * dim1_), dim0 (dim0_), dim1 (dim1_)
{
}

Matrix::Matrix (const unsigned & dim1_, const unsigned & dim0_, const double & value)
    : std::vector<double > (dim0_ * dim1_), dim0 (dim0_), dim1 (dim1_) 
{
  std::fill (this->begin(), this->end(), value);
}

Matrix::Matrix (const Matrix & mat)
    : std::vector<double > (mat), dim0 (mat.dim0), dim1 (mat.dim1) 
{
}

void Matrix::reinit (const unsigned & dim1_,
		     const unsigned & dim0_)
{
  dim0 = dim0_;
  dim1 = dim1_;
  this->resize (dim0 * dim1);
}

void Matrix::reinit (const unsigned & dim1_,
		     const unsigned & dim0_, 
		     const double &value)
{
  dim0 = dim0_;
  dim1 = dim1_;
  this->resize (dim0 * dim1);
  std::fill (this->begin(), this->end(), value);
}


// level 2
// v = beta * v + alpha * A^(t) * a
void MatrixOperation::mv (std::vector<double > & v,
			  const double & beta,
			  const double & alpha,
			  const Matrix & A, 
			  const bool & transA,
			  const std::vector<double > a)
{
  if (transA){
    assert (v.size() == A.d0());
    assert (a.size() == A.d1());
  }
  else{
    assert (v.size() == A.d1());
    assert (a.size() == A.d0());
  }
  char opA;
  transA ? (opA = 'n') : (opA = 't');
  // A is a matrix of dim m x n
  int m = A.d1();
  int n = A.d0();
  
  double al = alpha, be = beta;
  std::vector<double > Acpy (A);
  std::vector<double > acpy (a);
  int one = 1;
  
  dgemv (&opA, &n, &m, &al, &Acpy[0], &n, &acpy[0], &one, &be, &v[0], &one);
}



// v = alpha * A^(t) * a
void MatrixOperation::mv (std::vector<double > & v,
			  const double & alpha,
			  const Matrix & A, 
			  const bool & transA,
			  const std::vector<double > a)
{
  
  if (transA){
    v.resize (A.d0());
    assert (a.size() == A.d1());
  }
  else{
    v.resize (A.d1());
    assert (a.size() == A.d0());
  }
  std::fill (v.begin(), v.end(), 0.0);
  char opA;
  transA ? (opA = 'n') : (opA = 't');
  // A is a matrix of dim m x n
  int m = A.d1();
  int n = A.d0();
  
  double al = alpha, be = 0;
  std::vector<double > Acpy (A);
  std::vector<double > acpy (a);
  int one = 1;
  
  dgemv (&opA, &n, &m, &al, &Acpy[0], &n, &acpy[0], &one, &be, &v[0], &one);
}
  
// return u^(t) * A * v
double quadForm (const std::vector<double > & u,
		 const Matrix & A,
		 const std::vector<double > & v)
{
  std::vector<double > tmp;
  MatrixOperation::mv (tmp, 1., A, false, v);
  assert (u.size() == tmp.size());
  return std::inner_product (u.begin(), u.end(), tmp.begin(), 0.0);
}


// level 3
// C = alpha * A^(t) * B^(t)
void MatrixOperation::mm (Matrix & C,
			  const double & alpha,
			  const Matrix & A, const bool & transA,
			  const Matrix & B, const bool & transB)
{
  char opA, opB;
  transA ? (opA = 't') : (opA = 'n');
  transB ? (opB = 't') : (opB = 'n');
  
  int m, n, k;
  // A^transA is a matrix of dim m x k
  if (transA){
    m = A.d0();
    k = A.d1();
  }
  else {
    m = A.d1();
    k = A.d0();
  }
  // B^trans is a matrix of dim k x n
  if (transB){
    n = B.d1();
    assert (k == int(B.d0()));
  }
  else {
    n = B.d0();
    assert (k == int(B.d1()));
  }
  // as a result, C is a matrix of dim m x n
  C.reinit (m, n, 0.0);
  
  double al = alpha;
  double be = 0;
  std::vector<double > Acpy (A);
  std::vector<double > Bcpy (B);
  int ldA, ldB;
  transA ? (ldA = m) : (ldA = k);
  transB ? (ldB = k) : (ldB = n);

  dgemm (&opB, &opA, &n, &m, &k, &al, &Bcpy[0], &ldB, &Acpy[0], &ldA, &be, &C[0], &n);
} 


void MatrixOperation::print (const Matrix & A)
{
  std::vector<double >::const_iterator it = A.begin();
  for (unsigned i = 0; i < A.d1(); ++ i){
    for (unsigned j = 0; j < A.d0(); ++ j){
      std::cout << *(it++) << '\t';
    }
    std::cout << '\n';
  }
}

void MatrixOperation::print (const std::string & filename, const Matrix & A)
{
  FILE * fp = fopen (filename.c_str(), "w");
  std::vector<double >::const_iterator it = A.begin();
  for (unsigned i = 0; i < A.d1(); ++ i){
    for (unsigned j = 0; j < A.d0(); ++ j){
      fprintf(fp, "%.16e\t", *(it++));
    }
    fprintf (fp, "\n");
  }
  fclose(fp);
}
