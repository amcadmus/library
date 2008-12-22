#include "VectorOperation.h"

int main ()
{
  bool dim_valid;
  std::vector<double > A (9);
  std::vector<double > B (9);
  
  A[0] = 0;
  A[1] = 1;
  A[2] = 2;
  A[3] = 3;
  A[4] = 4;
  A[5] = 5;
  A[6] = 6;
  A[7] = 7;
  A[8] = 8;
  B[0] = 10 + 0;
  B[1] = 10 + 1;
  B[2] = 10 + 2;
  B[3] = 10 + 3;
  B[4] = 10 + 4;
  B[5] = 10 + 5;
  B[6] = 10 + 6;
  B[7] = 10 + 7;
  B[8] = 10 + 8;

  std::vector<double > a (3);
  a[0] = 9;
  a[1] = -3;
  a[2] = 45;
  
  std::vector<double > b (3);
  b[0] = -2.1;
  b[1] = 7.4;
  b[2] = 6.5;
  
  VectorOperation::print (a);
  VectorOperation::print (b);
  VectorOperation::printMatrix (A);

  std::cout << VectorOperation::quadForm(a, A, b) << std::endl;
  
  std::cout << VectorOperation::dot (a, a) << std::endl;
  
  
//   VectorOperation::mv (b, 5, 7, A, true, a);

//   VectorOperation::print (b);

//   std::vector<double > C (9, 1.);
//   VectorOperation::print (C, dim_valid);
//   VectorOperation::eyeCorrect (C, 1.);
//   VectorOperation::print (C, dim_valid);
  

//   VectorOperation::print (C, dim_valid);

//   VectorOperation::sadd (C, 3, 2, A, 5, B);
  
//   VectorOperation::mm (C, 3, A, true, B, false);
  
//   VectorOperation::print (A, dim_valid);
//   VectorOperation::print (B, dim_valid);
//   VectorOperation::print (C, dim_valid);

  return 0;
}
