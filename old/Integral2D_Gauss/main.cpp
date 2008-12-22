// #include "../Integral2D/Integral2D.h"
#include"Integral2D.h"

Integral2DMethod::Type Integral2DMethod::Simpson	= 0;
Integral2DMethod::Type Integral2DMethod::Gauss9		= 1;

static double global_hx;
static double global_hy;

template <class BinaryFunction, class ValueType >
ValueType Integral2D<BinaryFunction,ValueType >::cal_int (const Integral2DMethod::Type method ,
							  const BinaryFunction & f, 
							  const ValueType & xlo, const ValueType & xup,
							  const ValueType & ylo, const ValueType & yup,
							  const ValueType & tol, 
							  const int & nx,
							  const int & ny)
{
  hx = (xup - xlo) / ValueType (nx);
  hy = (yup - ylo) / ValueType (ny);
  global_hx = hx;
  global_hy = hy;
  
  int nboxes = nx * ny;
  ValueType nboxesi = 1./ nboxes;
  
  ValueType value = 0;
  for (int i = 0; i < nx; i ++){
    for (int j = 0; j < ny; j ++){
      ValueType x = i * hx;
      ValueType y = j * hy;
      box.reinit (f, x, x+hx, y, y+hy, tol * nboxesi);
      ValueType tmp;
      box.cal_int (tmp);

      value += tmp;
    }
  }
  return value;
}




template <class BinaryFunction, class ValueType >
void IntegralBox2D<BinaryFunction,ValueType >::reinit (const BinaryFunction & f_, 
						       const ValueType & xlo_, const ValueType & xup_,
						       const ValueType & ylo_, const ValueType & yup_,
						       const ValueType & tol_)
{
  f = f_;
  xlo = xlo_;
  ylo = ylo_;
  xup = xup_;
  yup = yup_;
  hx = (xup - xlo);
  hy = (yup - ylo);

  double sqrt3o5 = sqrt (0.6);
  double x0 = 0.5 * (1 + sqrt3o5) * xlo + 0.5 * (1 - sqrt3o5) * xup;
  double x1 = 0.5 * (xlo + xup);
  double x2 = 0.5 * (1 - sqrt3o5) * xlo + 0.5 * (1 + sqrt3o5) * xup;
  double y0 = 0.5 * (1 + sqrt3o5) * ylo + 0.5 * (1 - sqrt3o5) * yup;
  double y1 = 0.5 * (ylo + yup);
  double y2 = 0.5 * (1 - sqrt3o5) * ylo + 0.5 * (1 + sqrt3o5) * yup;
  
  coarse_value.resize (9);
  coarse_value[0] = f (x0, y0);
  coarse_value[1] = f (x1, y0);
  coarse_value[2] = f (x2, y0);
  coarse_value[3] = f (x0, y1);
  coarse_value[4] = f (x1, y1);
  coarse_value[5] = f (x2, y1);
  coarse_value[6] = f (x0, y2);
  coarse_value[7] = f (x1, y2);
  coarse_value[8] = f (x2, y2);
  
  coarse_int = 1./324. * (hx) * (hy) * (
      (coarse_value[0] + coarse_value[2] + coarse_value[6] + coarse_value[8]) * 25. +
      (coarse_value[1] + coarse_value[3] + coarse_value[5] + coarse_value[7]) * 40. +
      (coarse_value[4]) * 64. );
  
  tol = tol_;
}


template <class BinaryFunction, class ValueType >
void IntegralBox2D<BinaryFunction,ValueType >::cal_int (ValueType & value)
{
  ValueType hhx = 0.5 * hx;
  ValueType hhy = 0.5 * hy;
  
  IntegralBox2D box0;
  IntegralBox2D box1;
  IntegralBox2D box2;
  IntegralBox2D box3;

  box0.reinit (f, xlo, xlo+hhx, ylo, ylo+hhy, tol*0.25);
  box1.reinit (f, xlo+hhx, xup, ylo, ylo+hhy, tol*0.25);
  box2.reinit (f, xlo, xlo+hhx, ylo+hhy, yup, tol*0.25);
  box3.reinit (f, xlo+hhx, xup, ylo+hhy, yup, tol*0.25);
  
  value =  box0.coarseInt() + box1.coarseInt() + box2.coarseInt() + box3.coarseInt();

  ValueType diff = fabs (coarse_int - value);

  if (diff < tol * 255. ){ // stop
    return;
  }
  else if (hx < 1e-13 * global_hx || hy < 1e-13 * global_hy){
    return;
  }
  else {
    ValueType tmp;
    value = 0;
    box0.cal_int (tmp);
    value += tmp;
    box1.cal_int (tmp);
    value += tmp;
    box2.cal_int (tmp);
    value += tmp;
    box3.cal_int (tmp);
    value += tmp;
    
    return;
  }
}


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


typedef double (*F2) (double, double);
int main ()
{
  Integral2D<F2,double > integ;

  std::vector<double > R246 (15, 0.000);
  global_R246 = R246;
  
//   double value = integ.cal_int (f, 0, M_PI, 0, M_PI, 1e-6);
  double value = integ.cal_int (Integral2DMethod::Simpson, value0_12, 0, 1, 0, M_PI, 1e-20, 30, 40);
  
  std::cout.precision (16);
  std::cout << value << std::endl;
  
  return 0;
}

