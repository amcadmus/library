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

#include "Optimization.h"

void OPTIMIZATION_TARGET::
value (const std::vector<ScalorType > & x,
       ScalorType & value) const
{
  ScalorType sum = 0.;
  for (unsigned i = 0; i < x.size(); ++i){
    sum += x[i] * x[i];
  }
  value = 1-sum;
  value *= value;
}


void OPTIMIZATION_TARGET::
grad  (const std::vector<ScalorType > & x,
       std::vector<ScalorType > & g) const
{
  g.resize(x.size());
  for (unsigned i = 0; i < x.size(); ++i){
    std::vector<ScalorType > x1(x), x2(x);
    x1[i] += hg;
    x2[i] -= hg;
    ScalorType v1, v2;
    value(x1, v1);
    value(x2, v2);
    g[i] = (v1 - v2) / (2. * hg);
  }
}


void OPTIMIZATION_TARGET::
hessian (const std::vector<ScalorType > & x,
	 std::vector<ScalorType > & hessian) const
{
  unsigned stride = x.size();
  hessian.resize(stride * stride);
  for (unsigned i = 0; i < stride; ++i){
    for (unsigned j = 0; j < stride; ++j){
      std::vector<ScalorType > x1(x), x2(x), x3(x), x4(x);
      x1[i] += hh;
      x1[j] += hh;
      x2[i] += hh;
      x2[j] -= hh;
      x3[i] -= hh;
      x3[j] += hh;
      x4[i] -= hh;
      x4[j] -= hh;
      ScalorType v1, v2, v3, v4;
      value (x1, v1);
      value (x2, v2);
      value (x3, v3);
      value (x4, v4);
      hessian[i * stride + j] = (v1 - v2 - v3 + v4) / (4. * hh * hh);
    }
  }
}

int main(int argc, char * argv[])
{
  OPTIMIZATION_TARGET f;
  std::vector<ScalorType > x(3, 100.);
  x[2] = 50;
  x[0] = 200;
  
  OptimizationControler control;
  OptimizationCounter oc;
  
  SpeedestDescent sd;
  NewtonInteration ni;

  control.printDetail = true;
  sd.findMin (f, control, x, &oc);
  oc.print();
}
