#include "Matrix.h"
#include "../VectorOperation/VectorOperation.h"

void print (std::vector<double > & v )
{
  for (std::vector<double >::iterator pv = v.begin(); pv != v.end();){
    std::cout << *(pv ++) << std::endl;
  }
}

int main () 
{
  Matrix mat0 (3, 4);
  Matrix mat1 (3, 2);

  mat0 [0] = 1;
  mat0 [1] = 2;
  mat0 [2] = 3;
  mat0 [3] = 4;
  mat0 [4] = 5;
  mat0 [5] = 6;
  mat0 [6] = 7;
  mat0 [7] = 8;
  mat0 [8] = 9;
  mat0 [9] = 4;
  mat0 [10] = 5;
  mat0 [11] = 6;
  mat1 (0,0) = 9.5;
  mat1 (0,1) = 8.5;
  mat1 (1,0) = 7.5;
  mat1 (1,1) = 6.5;
  mat1 (2,0) = 5.5;
  mat1 (2,1) = 4.5;
  
  MatrixOperation::print (mat0);
  MatrixOperation::print (mat1);

  std::vector<double > a (3);
  a[0] = 9;
  a[1] = 8;
  a[2] = 7;
  
  std::vector<double > b (2);
  b[0] = 4.5;
  b[1] = 5.5;

//   MatrixOperation::mv (a, 1.5, 2.1, mat0, true, b);
  
//   print (a);

  Matrix C;
  MatrixOperation::mm (C, 2.2, mat0, true, mat1, false);

  MatrixOperation::print (C);
  

  
}

