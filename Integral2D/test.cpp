#include <iostream>
#include "Integral2D.h"
#include <cmath>

typedef double (*F2) (double , double);

double indx = 50;
std::vector<double > global_R246;

double f (double x, double y)
{
  return sin(x) * sin(y);
}

double value11_1 (double r, double theta)
{
  double COS = cos(theta);
  double SIN = sin(theta);
  double COS2 = COS * COS;
  double COS4 = COS2 * COS2;
  double SIN2 = SIN * SIN;
  double SIN4 = SIN2 * SIN2;
  double r2 = r * r;
  double r4 = r2 * r2;
  double R246Q2 = r2 * (global_R246[0] * COS2 + global_R246[1] * COS * SIN + global_R246[2] * SIN2);
  double R246Q4 = r2 * r2 * (global_R246[0+3] * COS4 + global_R246[1+3] * COS2 * COS * SIN + global_R246[2+3] * COS2 * SIN2 + global_R246[3+3] * COS * SIN * SIN2 + global_R246[4+3] * SIN4);
  double R246Q6 = r2 * r2 * r2 * (global_R246[0+8] * COS2 * COS2 * COS2 + global_R246[1+8] * COS2 * COS2 * COS * SIN + global_R246[2+8] * COS2 * COS2 * SIN2 + global_R246[3+8] * COS2 * COS * SIN * SIN2 + global_R246[4+8] * COS2 * SIN2 * SIN2 + global_R246[5+8] * COS * SIN * SIN2 * SIN2 + global_R246[6+8] * SIN2 * SIN2 * SIN2);
  return r4 * r4 * r4 * COS4 * COS4 * COS2 * COS * SIN * r * pow(1-r2,indx) * exp (R246Q2 + R246Q4 + R246Q6);
}


double value0_12 (double r, double theta)
{
  double COS = cos(theta);
  double SIN = sin(theta);
  double COS2 = COS * COS;
  double COS4 = COS2 * COS2;
  double SIN2 = SIN * SIN;
  double SIN4 = SIN2 * SIN2;
  double r2 = r * r;
  double r4 = r2 * r2;
  double R246Q2 = r2 * (global_R246[0] * COS2 + global_R246[1] * COS * SIN + global_R246[2] * SIN2);
  double R246Q4 = r2 * r2 * (global_R246[0+3] * COS4 + global_R246[1+3] * COS2 * COS * SIN + global_R246[2+3] * COS2 * SIN2 + global_R246[3+3] * COS * SIN * SIN2 + global_R246[4+3] * SIN4);
  double R246Q6 = r2 * r2 * r2 * (global_R246[0+8] * COS2 * COS2 * COS2 + global_R246[1+8] * COS2 * COS2 * COS * SIN + global_R246[2+8] * COS2 * COS2 * SIN2 + global_R246[3+8] * COS2 * COS * SIN * SIN2 + global_R246[4+8] * COS2 * SIN2 * SIN2 + global_R246[5+8] * COS * SIN * SIN2 * SIN2 + global_R246[6+8] * SIN2 * SIN2 * SIN2);
  return r4 * r4 * r4 * SIN4 * SIN4 * SIN4 * r * pow(1-r2,indx) * exp (R246Q2 + R246Q4 + R246Q6);
}

double value20 (double r, double theta)
{
  double COS = cos(theta);
  double SIN = sin(theta);
  double COS2 = COS * COS;
  double SIN2 = SIN * SIN;
  double r2 = r * r;
  double R246Q2 = r2 * (global_R246[0] * COS2 + global_R246[1] * COS * SIN + global_R246[2] * SIN2);
  double R246Q4 = r2 * r2 * (global_R246[0+3] * COS2 * COS2 + global_R246[1+3] * COS2 * COS * SIN + global_R246[2+3] * COS2 * SIN2 + global_R246[3+3] * COS * SIN * SIN2 + global_R246[4+3] * SIN2 * SIN2);
  double R246Q6 = r2 * r2 * r2 * (global_R246[0+8] * COS2 * COS2 * COS2 + global_R246[1+8] * COS2 * COS2 * COS * SIN + global_R246[2+8] * COS2 * COS2 * SIN2 + global_R246[3+8] * COS2 * COS * SIN * SIN2 + global_R246[4+8] * COS2 * SIN2 * SIN2 + global_R246[5+8] * COS * SIN * SIN2 * SIN2 + global_R246[6+8] * SIN2 * SIN2 * SIN2);
  return r2 * COS2 * r * pow(1-r2,indx) * exp (R246Q2 + R246Q4 + R246Q6);
}


int main ()
{
  Integral2D<F2,double > integ;

  std::vector<double > R246 (15, 0.000);
  global_R246 = R246;
  
  double value0;
  double value1;
  double e = 2.9615315361867009e-04;
  // 1.914677683422673e-10;
//   e = 4;
  double tol = 1e-14;
  
  for (int i = 0; i < 1; i ++){
//     value0 = integ.cal_int (Integral2DInfo::Gauss9, f, 0, M_PI, 0, M_PI, 1e-10);
//     value1 = integ.cal_int (Integral2DInfo::Simpson, f, 0, M_PI, 0, M_PI, 1e-8);
    value0 = integ.cal_int (Integral2DInfo::Gauss9, value20, 0, 1, 0, M_PI, tol, 3, 4);
    value1 = integ.cal_int (Integral2DInfo::Simpson, value20, 0, 1, 0, M_PI, tol, 3, 4);
  }
  
  
  std::cout.precision (16);
  std::cout << std::scientific << value0 - e << std::endl;
  std::cout << std::scientific << value1 - e << std::endl;
//   std::cout << integ.cal_int (Integral2DInfo::Gauss9, value20, 0, 1, 0, M_PI, 1e-17, 3, 4) << std::endl;
  
  
  return 0;
}
