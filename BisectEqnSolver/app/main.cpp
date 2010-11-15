#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>

#include "BisectEqnSolver.h"

// typedef  double VARIABLETYPE;
// typedef double VALUETYPE;
typedef double (*TYPE_PDF) (double);


double f (double x)
{
  return sin(x);
}


int main(int argc, char * argv[])
{
  BisectEqnSolver<double, double, TYPE_PDF> solver;
  double x;
  solver.solve (f, 2., 4., 1e-7, x);
  std::cout << "result is " << x - M_PI << std::endl;

  return 0;
}
