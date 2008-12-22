#include "VectorOperation.h"

void VectorOperation::zero (std::vector<double > & v)
{
  std::fill (v.begin(), v.end(), 0.);
}

void VectorOperation::scale (std::vector<double > & v, const double & scale)
{
  std::vector<double >::iterator pv = v.begin();
  for (; pv != v.end(); ++ pv){
    *pv *= scale;
  }
}

void VectorOperation::add (std::vector<double > & v, 
			   const double & alpha, const std::vector<double > & a)
{
  assert (v.size() == a.size());
  std::vector<double >::iterator pv = v.begin();
  std::vector<double >::const_iterator pa = a.begin();
  for (; pv != v.end(); ++ pv, ++ pa){
    *pv += alpha * *pa;
  }
}

void VectorOperation::add (std::vector<double > & v,
			   const double & alpha, const std::vector<double > & a,
			   const double & beta,  const std::vector<double > & b)
{
  assert (v.size() == a.size());
  assert (v.size() == b.size());
  std::vector<double >::iterator pv = v.begin();
  std::vector<double >::const_iterator pa = a.begin();
  std::vector<double >::const_iterator pb = b.begin();
  for (; pv != v.end(); ++pv, ++pa, ++pb){
    *pv += alpha * *pa + beta * *pb;
  }
}

void VectorOperation::sadd (std::vector<double > & v, const double & s, 
			    const double & alpha, const std::vector<double > & a)
{
  assert (v.size() == a.size());
  std::vector<double >::iterator pv = v.begin();
  std::vector<double >::const_iterator pa = a.begin();
  for (; pv != v.end(); ++ pv, ++ pa){
    *pv = s * *pv + alpha * *pa;
  }
}
 
void VectorOperation::sadd (std::vector<double > & v, const double & s,
			    const double & alpha, const std::vector<double > & a,
			    const double & beta,  const std::vector<double > & b)
{
  assert (v.size() == a.size());
  assert (v.size() == b.size());
  std::vector<double >::iterator pv = v.begin();
  std::vector<double >::const_iterator pa = a.begin();
  std::vector<double >::const_iterator pb = b.begin();
  for (; pv != v.end(); ++pv, ++pa, ++pb){
    *pv = s * *pv + alpha * *pa + beta * *pb;
  }
}

double VectorOperation::dot (const std::vector<double > & u, const std::vector<double > & v)
{
  assert (v.size() == u.size());
  return std::inner_product (u.begin(), u.end(), v.begin(), 0.);
}

unsigned VectorOperation::dim (const std::vector<double > & A, bool & dim_valid)
{
  double sqrtAsize = sqrt (A.size());
  unsigned n = unsigned (sqrtAsize);
  dim_valid = (double(n) == sqrtAsize);
  return n;
}

void VectorOperation::print (const std::vector<double > & v, const bool & inColume)
{
  if (inColume){
    for (std::vector<double >::const_iterator pv = v.begin(); pv != v.end();){
      std::cout << *(pv ++) << std::endl;
    }
  }
  else{
    for (std::vector<double >::const_iterator pv = v.begin(); pv != v.end();){
      std::cout << *(pv ++) << "\t";
    }
    std::cout << std::endl;
  }
}

void VectorOperation::print (const std::string & filename, const std::vector<double > & v)
{
  FILE * fp = fopen (filename.c_str(), "w");
  std::vector<double >::const_iterator it = v.begin();
  for (; it != v.end();){
    fprintf(fp, "%.16e\n", *(it++));
  }
  fclose (fp);
} 

bool VectorOperation::read (const std::string & filename, std::vector<double > & v)
{
  FILE * fp = fopen (filename.c_str(), "r");
  if (fp == NULL) return false;
  double tmp;
  v.clear ();
  int info;
  
  while ((info = fscanf (fp, "%lf", &tmp)) != EOF){
    assert (info == 1);
    v.push_back (tmp);
  }
  
  fclose (fp);
  return true;
}

void VectorOperation::mv (std::vector<double > & v,
			  const double & beta,
			  const double & alpha, const std::vector<double > & A, const bool & transA, const std::vector<double > & a)
{
  char opA;
  transA ? (opA = 'n') : (opA = 't');
  bool dim_valid;
  unsigned n = dim (A, dim_valid);
  assert (dim_valid);
  assert (n == a.size());
  assert (n == v.size());
  int m = int(n);

  double al = alpha;
  double be = beta;
  std::vector<double > Acpy (A);
  std::vector<double > acpy (a);
  int one = 1;
  dgemv (&opA, &m, &m, &al, &Acpy[0], &m, &acpy[0], &one, &be, &v[0], &one);
}

void VectorOperation::mv (std::vector<double > & v,
			  const double & alpha, const std::vector<double > & A, const bool & transA, const std::vector<double > & a)
{
  char opA;
  transA ? (opA = 'n') : (opA = 't');
  bool dim_valid;
  unsigned nn = dim (A, dim_valid);
  assert (dim_valid);
  assert (nn == a.size());
  v.resize (nn);
  std::fill (v.begin(), v.end(), 0.0);
  int n = nn;

  double al = alpha;
  double be = 0.;
  std::vector<double > Acpy (A);
  std::vector<double > acpy (a);
  int one = 1;
  dgemv (&opA, &n, &n, &al, &Acpy[0], &n, &acpy[0], &one, &be, &v[0], &one);
}

double VectorOperation::quadForm (const std::vector<double > & u, const std::vector<double > & A, std::vector<double > & v)
{
  std::vector<double > tmp;
  mv (tmp, 1., A, false, v);
  return dot (tmp, u);
}

void VectorOperation::eyeCorrect (std::vector<double > & A, const double & alpha)
{
  bool dim_valid;
  int n = dim (A, dim_valid);
  assert (dim_valid);

  std::vector<double >::iterator pa = A.begin();
  for (int i = 0; i < n; i ++){
    *pa += alpha;
    pa += n + 1;
  }
}


void VectorOperation::mm (std::vector<double > & C,
			  const double & alpha, 
			  const std::vector<double > & A, const bool & transA, 
			  const std::vector<double > & B, const bool & transB)
{
  char opA, opB;
  transA ? (opA = 't') : (opA = 'n');
  transB ? (opB = 't') : (opB = 'n');
  
  double sqrtAsize = sqrt (A.size());
  int n = int (sqrtAsize);
  assert (n == sqrtAsize);
  assert (n * n == int(B.size()));
  C.resize (n * n);
  std::fill (C.begin(), C.end(), 0.0);
  
  double al = alpha;
  double be = 0;
  std::vector<double > Acpy (A);
  std::vector<double > Bcpy (B);
  dgemm (&opB, &opA, &n, &n, &n, &al, &Bcpy[0], &n, &Acpy[0], &n, &be, &C[0], &n);
}

void VectorOperation::printMatrix (const std::vector<double > & A)
{
  bool valid;
  int n = dim (A, valid);
  assert (valid);
  
  std::vector<double >::const_iterator it = A.begin();
  for (int i = 0; i < n; ++ i){
    for (int j = 0; j < n; ++ j){
      std::cout << *(it ++) << '\t';
    }
    std::cout << '\n';
  }
}
 

void VectorOperation::printMatrix (const std::string & filename, 
				   const std::vector<double > & A)
{
  bool valid;
  int n = dim (A, valid);
  assert (valid);
 
  FILE * fp = fopen (filename.c_str(), "w");
  std::vector<double >::const_iterator it = A.begin();
  for (int i = 0; i < n; ++ i){
    for (int j = 0; j < n; ++ j){
      fprintf(fp, "%.16e\t", *(it++));
    }
    fprintf (fp, "\n");
  }
  fclose(fp);
}

void VectorOperation::printMatrix (const std::vector<double > & A,
				   const int & m, 
				   const int & n )
{
  std::vector<double >::const_iterator it = A.begin();
  for (int i = 0; i < m; ++ i){
    for (int j = 0; j < n; ++ j){
      std::cout << *(it ++) << '\t';
    }
    std::cout << ";\n";
  }
}

void VectorOperation::printMatrix (const std::string & filename, 
				   const std::vector<double > & A,
				   const int & m, 
				   const int & n )
{
  FILE * fp = fopen (filename.c_str(), "w");
  std::vector<double >::const_iterator it = A.begin();
  for (int i = 0; i < m; ++ i){
    for (int j = 0; j < n; ++ j){
      fprintf(fp, "%.16e\t", *(it++));
    }
    fprintf (fp, "\n");
  }
  fclose(fp);
}

void VectorOperation::printVector (const std::vector<double > & v)
{
  for (std::vector<double >::const_iterator pv = v.begin(); 
       pv != v.end(); ++pv){
    printf ("%.16e\n", *pv);
  }
}

void VectorOperation::printVector (const std::vector<double > & u,
				   const std::vector<double > & v)
{
  std::vector<double >::const_iterator pv = v.begin();
  std::vector<double >::const_iterator pu = u.begin();
  for ( ;pu != u.end(); ++pv, ++pu){
    printf ("%.16e\t%.16e\n", *pu, *pv);
  }
}

void VectorOperation::printVector (const std::string & filename, 
				   const std::vector<double > & v)
{
  FILE * fp = fopen (filename.c_str(), "w");
  for (std::vector<double >::const_iterator pv = v.begin(); 
       pv != v.end(); ++pv){
    fprintf (fp, "%.16e\n", *pv);
  }
  fclose(fp);
}

void VectorOperation::printVector (const std::string & filename, 
				   const std::vector<double > & u,
				   const std::vector<double > & v)
{
  FILE * fp = fopen (filename.c_str(), "w");
  std::vector<double >::const_iterator pv = v.begin();
  std::vector<double >::const_iterator pu = u.begin();
  for ( ;pu != u.end(); ++pv, ++pu){
    fprintf (fp, "%.16e\t%.16e\n", *pu, *pv);
  }
  fclose(fp);
}

void VectorOperation::loadVector (const std::string & filename, 
				  std::vector<double > & v )
{
  FILE * fp = fopen (filename.c_str(), "r");
  if (fp == NULL ){
    std::cerr << "cannot open file " << filename << std::endl;
    return ;
  }
  v.clear();
  double tmp;
  while (fscanf (fp, "%lf", &tmp) != EOF){
    v.push_back (tmp);
  }
  fclose (fp);
}

void VectorOperation::loadVector (const std::string & filename, 
				  std::vector<double > & u ,
				  std::vector<double > & v )
{
  FILE * fp = fopen (filename.c_str(), "r");
  if (fp == NULL ){
    std::cerr << "cannot open file " << filename << std::endl;
    return ;
  }
  v.clear();
  u.clear();
  double tmp0;
  double tmp1;
  while (fscanf (fp, "%lf%lf", &tmp0, &tmp1) != EOF){
    u.push_back (tmp0);
    v.push_back (tmp1);
  }
  fclose (fp);
}
