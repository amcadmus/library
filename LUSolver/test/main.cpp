#include"LUSolver.h"

int main ()
{
  std::vector<double > matrix (9);
  matrix[0] = 10;
  matrix[1] = 5;
  matrix[2] = 2;
  matrix[3] = 5;
  matrix[4] = 20;
  matrix[5] = 4;
  matrix[6] = 2;
  matrix[7] = 4;
  matrix[8] = 15;

  std::vector<double > rhs (3);
  rhs[0] = -1;
  rhs[1] = 5;
  rhs[2] = -10;
  std::vector<double > result;

  LUSolver<double > solver;
  solver.reinit (matrix, 3);
  solver.solve (rhs, result);
  
  std::cout.precision (16);
  std::cout << std::scientific 
	    << result[0] << "\n"
	    << result[1] << "\n"
	    << result[2] << std::endl;
  
  return 0;
}

  
