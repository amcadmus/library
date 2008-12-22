#include "LQFactorization.h"

int LQFactorization<double >::run ( const std::vector<double > & matrix, 
				    const int & mm, // mm < nn
				    const int & nn,
				    std::vector<double > & L,
				    std::vector<double > & Q)
{
  Q.clear();
  Q.resize (nn * nn);
  std::copy (matrix.begin(), matrix.end(), Q.begin());
//   std::vector<double > workmat (matrix);
  std::vector<double > work (1);
  int lwork = -1;
  int m = nn;
  int n = mm;
  int lda = m;
  std::vector<double > tau (m);
  int info;
  
  dgeqrf (&m, &n, &Q[0], &lda, &tau[0], &work[0], &lwork, &info);
  lwork = 64;
  while (lwork < work[0]) lwork += 64;
  work.clear ();
  work.resize (lwork);
//  std::cout << "lwork is " << lwork << "\nwork[0] is " << work[0] << std::endl;
  
  dgeqrf (&m, &n, &Q[0], &lda, &tau[0], &work[0], &lwork, &info);
  if (info != 0) return info;
  L.clear();
  L.resize (n*n, 0);
  std::vector<double >::iterator pL = L.begin();
  std::vector<double >::iterator pQ = Q.begin();
  for (int i = 0; i < n; ++i ){
    for (int j = 0; j < i+1 ; ++j)
      *(pL++) = *(pQ++);
    pL += n-i-1;
    pQ += m-i-1;
  }
  
  lwork = -1;
  dorgqr (&m, &m, &n, &Q[0], &lda, &tau[0], &work[0], &lwork, &info);
  lwork = 64;
  while (lwork < work[0]) lwork += 64;
  work.clear ();
  work.resize (lwork);
//  std::cout << "lwork is " << lwork << "\nwork[0] is " << work[0] << std::endl;
  
  dorgqr (&m, &m, &n, &Q[0], &lda, &tau[0], &work[0], &lwork, &info);
  // if (info != 0) return info;
  return info;
}

// void printMatrix (const std::vector<double > & A,
// 		  int m, int n)
// {
//   std::vector<double >::const_iterator it = A.begin();
//   for (int i = 0; i < m; ++ i){
//     for (int j = 0; j < n; ++ j){
//       std::cout << *(it ++) << '\t';
//     }
//     std::cout << ";\n";
//   }
// }


// int main(int argc, char * argv[])
// {
//   LQFactorization<double > lq ;
// //   double a[9] = {1, 2, 5, 6, 3, 4, 7, 8};
// //   int m = 2;
// //   int n = 4;
//   double a[13] = {1,2,3,10,4,5,6,11,7,8,9,12};
//   int m = 3;
//   int n = 4;
//   std::vector<double > aa(m*n);
//   std::copy (&a[0], &a[m*n], aa.begin());
  
//   std::vector<double > L, Q;
//   int info = lq.run (aa, m, n, L, Q);
//   std::cout << info << std::endl ;
//   std::cout << std::endl;
  
//   std::cout.precision (16);
//   printMatrix (aa, m, n);
//   std::cout << std::endl;
//   printMatrix (L, m, m);
//   std::cout << std::endl;
//   printMatrix (Q, n, n);
//   std::cout << std::endl;
  
// }
