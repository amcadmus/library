#include "Integral2D.h"
#include <cmath>

Integral2DInfo::Method Integral2DInfo::Simpson	= 0;
Integral2DInfo::Method Integral2DInfo::Gauss9	= 1;
Integral2DInfo::Method Integral2DInfo::Gauss16	= 2;

static double global_hx;
static double global_hy;

template <class BinaryFunction, class ValueType >
ValueType Integral2D<BinaryFunction,ValueType >::cal_int (const Integral2DInfo::Method & method,
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

  if (method == Integral2DInfo::Simpson){
    for (int i = 0; i < nx; i ++){
      for (int j = 0; j < ny; j ++){
	ValueType x = i * hx + xlo;
	ValueType y = j * hy + ylo;
	box_Simpson.reinit (f, x, x+hx, y, y+hy, tol * nboxesi);
	ValueType tmp;
	box_Simpson.cal_int (tmp);

	value += tmp;
      }
    }
  }
  
  else if (method == Integral2DInfo::Gauss9){
    for (int i = 0; i < nx; i ++){
      for (int j = 0; j < ny; j ++){
	ValueType x = i * hx + xlo;
	ValueType y = j * hy + ylo;
	box_Gauss9.reinit (f, x, x+hx, y, y+hy, tol * nboxesi);
	ValueType tmp;
	box_Gauss9.cal_int (tmp);

	value += tmp;
      }
    }
  }
  else if (method == Integral2DInfo::Gauss16){
    for (int i = 0; i < nx; i ++){
      for (int j = 0; j < ny; j ++){
	ValueType x = i * hx + xlo;
	ValueType y = j * hy + ylo;
	box_Gauss16.reinit (f, x, x+hx, y, y+hy, tol * nboxesi);
	ValueType tmp;
	box_Gauss16.cal_int (tmp);

	value += tmp;
      }
    }
  }
  else{
    std::cerr << "Wrong integration method\n";
  }
  
  return value;
}


////////////////////////////////////////////////////////////////////////////////////////////////////


template <class BinaryFunction, class ValueType >
void IntegralBox2DSimpson<BinaryFunction,ValueType >::reinit (const BinaryFunction & f_, 
							      const ValueType & xlo_, const ValueType & xup_,
							      const ValueType & ylo_, const ValueType & yup_,
							      const ValueType & tol_)
{
  f = f_;
  xlo = xlo_;
  ylo = ylo_;
  xup = xup_;
  yup = yup_;
  hx = 0.5 * (xup - xlo);
  hy = 0.5 * (yup - ylo);
  one_36 = 1./36.;

  coarse_value.resize (9);
  coarse_value[0] = f (xlo,	ylo	);
  coarse_value[1] = f (xlo + hx,ylo	);
  coarse_value[2] = f (xup,	ylo	);
  coarse_value[3] = f (xlo,	ylo + hy);
  coarse_value[4] = f (xlo + hx,ylo + hy);
  coarse_value[5] = f (xup,	ylo + hy);
  coarse_value[6] = f (xlo,	yup	);
  coarse_value[7] = f (xlo + hx,yup	);
  coarse_value[8] = f (xup,	yup	);
  
  coarse_int = ( one_36 * (2.*hx) * (2.*hy) *
		 ((coarse_value[0] + coarse_value[2] + coarse_value[6] + coarse_value[8]) +
		  (coarse_value[1] + coarse_value[3] + coarse_value[5] + coarse_value[7]) * 4. +
		  (coarse_value[4]) * 16. )
      );
  tol = tol_;
}

template <class BinaryFunction, class ValueType >
void IntegralBox2DSimpson<BinaryFunction,ValueType >::reinit (const BinaryFunction & f_, 
							      const ValueType & xlo_, const ValueType & xup_, 
							      const ValueType & ylo_, const ValueType & yup_,
							      const std::vector<ValueType > & coarse_value_,
							      const ValueType & coarse_int_,
							      const ValueType & tol_)
{
  f = f_;
  xlo = xlo_;
  ylo = ylo_;
  xup = xup_;
  yup = yup_;
  hx = 0.5 * (xup - xlo);
  hy = 0.5 * (yup - ylo);
  one_36 = 1./36.;

  coarse_int = coarse_int_;
  coarse_value = coarse_value_;
  
  tol = tol_;
}

template <class BinaryFunction, class ValueType >
void IntegralBox2DSimpson<BinaryFunction,ValueType >::cal_int (ValueType & value)
{
  ValueType hhx = 0.5 * hx;
  ValueType hhy = 0.5 * hy;
//   if ( (hhx < 1e-13 || hhy < 1e-13) && hhx * hhy * coarse_value[4] < 1e-30 ){
// //   if (hhx < 1e-13 || hhy < 1e-13){
// //     std::cout << hhx << " " <<  coarse_value[4] <<  " Over refine! return!\n" ;
//     value = coarse_int;
//     return;
//   }
  
  std::vector<ValueType > refine_value0 (9);
  std::vector<ValueType > refine_value1 (9);
  std::vector<ValueType > refine_value2 (9);
  std::vector<ValueType > refine_value3 (9);
  
  ValueType refine_int0;
  ValueType refine_int1;
  ValueType refine_int2;
  ValueType refine_int3;
  
  refine_value0[0] = coarse_value[0];
  refine_value0[2] = coarse_value[1];
  refine_value0[6] = coarse_value[3];
  refine_value0[8] = coarse_value[4];
  refine_value0[1] = f (xlo+hhx,	ylo		);
  refine_value0[3] = f (xlo,		ylo+hhy		);
  refine_value0[4] = f (xlo+hhx,	ylo+hhy		);
  refine_value0[5] = f (xlo+hhx*2.,	ylo+hhy		);
  refine_value0[7] = f (xlo+hhx,	ylo+hhy*2.	);
  
  refine_int0 = one_36 * ( (refine_value0[0] + refine_value0[2] + refine_value0[6] + refine_value0[8]) +
			   (refine_value0[1] + refine_value0[3] + refine_value0[5] + refine_value0[7]) * 4. +
			   (refine_value0[4]) * 16. ) * hx * hy;
  

  refine_value1[0] = coarse_value[1];
  refine_value1[2] = coarse_value[2];
  refine_value1[6] = coarse_value[4];
  refine_value1[8] = coarse_value[5];
  refine_value1[3] = refine_value0[5];
  refine_value1[1] = f (xlo+hhx*3.,	ylo		);
  refine_value1[4] = f (xlo+hhx*3.,	ylo+hhy		);
  refine_value1[5] = f (xup,		ylo+hhy		);
  refine_value1[7] = f (xlo+hhx*3.,	ylo+hhy*2.	);
  
  refine_int1 = one_36 * ( (refine_value1[0] + refine_value1[2] + refine_value1[6] + refine_value1[8]) +
			   (refine_value1[1] + refine_value1[3] + refine_value1[5] + refine_value1[7]) * 4. +
			   (refine_value1[4]) * 16. ) * hx * hy;


  refine_value2[0] = coarse_value[3];
  refine_value2[2] = coarse_value[4];
  refine_value2[6] = coarse_value[6];
  refine_value2[8] = coarse_value[7];
  refine_value2[1] = refine_value0[7];
  refine_value2[3] = f (xlo,		ylo+hhy*3.	);
  refine_value2[4] = f (xlo+hhx,	ylo+hhy*3.	);
  refine_value2[5] = f (xlo+hhx*2.,	ylo+hhy*3.	);
  refine_value2[7] = f (xlo+hhx,	yup		);

  refine_int2 = one_36 * ( (refine_value2[0] + refine_value2[2] + refine_value2[6] + refine_value2[8]) +
			   (refine_value2[1] + refine_value2[3] + refine_value2[5] + refine_value2[7]) * 4. +
			   (refine_value2[4]) * 16. ) * hx * hy;

  refine_value3[0] = coarse_value[4];
  refine_value3[2] = coarse_value[5];
  refine_value3[6] = coarse_value[7];
  refine_value3[8] = coarse_value[8];
  refine_value3[1] = refine_value1[7];
  refine_value3[3] = refine_value2[5];
  refine_value3[4] = f (xlo+hhx*3.,	ylo+hhy*3.	);
  refine_value3[5] = f (xup,		ylo+hhy*3.	);
  refine_value3[7] = f (xlo+hhx*3.,	yup		);
  
  refine_int3 = one_36 * ( (refine_value3[0] + refine_value3[2] + refine_value3[6] + refine_value3[8]) +
			   (refine_value3[1] + refine_value3[3] + refine_value3[5] + refine_value3[7]) * 4. +
			   (refine_value3[4]) * 16. ) * hx * hy;

  value =  refine_int0 + refine_int1 + refine_int2 + refine_int3;
  ValueType diff = fabs (coarse_int - value);

  if (diff < tol * 63. ){ // stop
    value = (value - 1./16.*coarse_int) / (1 - 1./16.);
    return ;
  }
  else if (hx < 1e-13 * global_hx || hy < 1e-13 * global_hy){
    return ;
  }
  else{
    IntegralBox2DSimpson box0;
    IntegralBox2DSimpson box1;
    IntegralBox2DSimpson box2;
    IntegralBox2DSimpson box3;
    
    box0.reinit (f, xlo, xlo+hx, ylo, ylo+hy, refine_value0, refine_int0, tol*0.25);
    box1.reinit (f, xlo+hx, xup, ylo, ylo+hy, refine_value1, refine_int1, tol*0.25);
    box2.reinit (f, xlo, xlo+hx, ylo+hy, yup, refine_value2, refine_int2, tol*0.25);
    box3.reinit (f, xlo+hx, xup, ylo+hy, yup, refine_value3, refine_int3, tol*0.25);
    
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


template <class BinaryFunction, class ValueType >
void IntegralBox2DGauss9<BinaryFunction,ValueType >::reinit (const BinaryFunction & f_, 
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

  ValueType sqrt3o5 = sqrt (0.6);
  ValueType x0 = 0.5 * (1 + sqrt3o5) * xlo + 0.5 * (1 - sqrt3o5) * xup;
  ValueType x1 = 0.5 * (xlo + xup);
  ValueType x2 = 0.5 * (1 - sqrt3o5) * xlo + 0.5 * (1 + sqrt3o5) * xup;
  ValueType y0 = 0.5 * (1 + sqrt3o5) * ylo + 0.5 * (1 - sqrt3o5) * yup;
  ValueType y1 = 0.5 * (ylo + yup);
  ValueType y2 = 0.5 * (1 - sqrt3o5) * ylo + 0.5 * (1 + sqrt3o5) * yup;
  
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
void IntegralBox2DGauss9<BinaryFunction,ValueType >::cal_int (ValueType & value)
{
  ValueType hhx = 0.5 * hx;
  ValueType hhy = 0.5 * hy;
  
  IntegralBox2DGauss9 box0;
  IntegralBox2DGauss9 box1;
  IntegralBox2DGauss9 box2;
  IntegralBox2DGauss9 box3;

  box0.reinit (f, xlo, xlo+hhx, ylo, ylo+hhy, tol*0.25);
  box1.reinit (f, xlo+hhx, xup, ylo, ylo+hhy, tol*0.25);
  box2.reinit (f, xlo, xlo+hhx, ylo+hhy, yup, tol*0.25);
  box3.reinit (f, xlo+hhx, xup, ylo+hhy, yup, tol*0.25);
  
  value =  box0.coarseInt() + box1.coarseInt() + box2.coarseInt() + box3.coarseInt();

  ValueType diff = fabs (coarse_int - value);

  if (diff < tol * 255. ){ // stop
//     value = (value - coarse_int/64.) / (1 - 1./64.);
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



template <class BinaryFunction, class ValueType >
void IntegralBox2DGauss16<BinaryFunction,ValueType >::reinit (const BinaryFunction & f_, 
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

  ValueType value0 = sqrt ((3 - 2 * sqrt(1.2)) / 7.);
  ValueType value1 = sqrt ((3 + 2 * sqrt(1.2)) / 7.);
  
  ValueType x0 = 0.5 * (xup - xlo) * (-value1) + 0.5 * (xup + xlo);
  ValueType x1 = 0.5 * (xup - xlo) * (-value0) + 0.5 * (xup + xlo);
  ValueType x2 = 0.5 * (xup - xlo) * (+value0) + 0.5 * (xup + xlo);
  ValueType x3 = 0.5 * (xup - xlo) * (+value1) + 0.5 * (xup + xlo);
  ValueType y0 = 0.5 * (yup - ylo) * (-value1) + 0.5 * (yup + ylo);
  ValueType y1 = 0.5 * (yup - ylo) * (-value0) + 0.5 * (yup + ylo);
  ValueType y2 = 0.5 * (yup - ylo) * (+value0) + 0.5 * (yup + ylo);
  ValueType y3 = 0.5 * (yup - ylo) * (+value1) + 0.5 * (yup + ylo);
  
  coarse_value.resize (16);
  coarse_value[0] = f (x0, y0);
  coarse_value[1] = f (x1, y0);
  coarse_value[2] = f (x2, y0);
  coarse_value[3] = f (x3, y0);
  coarse_value[4] = f (x0, y1);
  coarse_value[5] = f (x1, y1);
  coarse_value[6] = f (x2, y1);
  coarse_value[7] = f (x3, y1);
  coarse_value[8] = f (x0, y2);
  coarse_value[9] = f (x1, y2);
  coarse_value[10] = f (x2, y2);
  coarse_value[11] = f (x3, y2);
  coarse_value[12] = f (x0, y3);
  coarse_value[13] = f (x1, y3);
  coarse_value[14] = f (x2, y3);
  coarse_value[15] = f (x3, y3);
  
  double weigh0 = (18 + sqrt(30.)) / 36.;
  double weigh1 = (18 - sqrt(30.)) / 36.;
  
  coarse_int = 1./4. * (hx) * (hy) * (
      weigh1 * weigh1 * (coarse_value[0] + coarse_value[3] + coarse_value[12] + coarse_value[15]) +
      weigh0 * weigh1 * (coarse_value[1] + coarse_value[2] + coarse_value[13] + coarse_value[14] + coarse_value[4] + coarse_value[8] + coarse_value[7] + coarse_value[11]) +
      weigh0 * weigh0 * (coarse_value[5] + coarse_value[6] + coarse_value[9] + coarse_value[10]) );
  
  tol = tol_;
}


template <class BinaryFunction, class ValueType >
void IntegralBox2DGauss16<BinaryFunction,ValueType >::cal_int (ValueType & value)
{
  ValueType hhx = 0.5 * hx;
  ValueType hhy = 0.5 * hy;
  
  IntegralBox2DGauss16 box0;
  IntegralBox2DGauss16 box1;
  IntegralBox2DGauss16 box2;
  IntegralBox2DGauss16 box3;

  box0.reinit (f, xlo, xlo+hhx, ylo, ylo+hhy, tol*0.25);
  box1.reinit (f, xlo+hhx, xup, ylo, ylo+hhy, tol*0.25);
  box2.reinit (f, xlo, xlo+hhx, ylo+hhy, yup, tol*0.25);
  box3.reinit (f, xlo+hhx, xup, ylo+hhy, yup, tol*0.25);
  
  value =  box0.coarseInt() + box1.coarseInt() + box2.coarseInt() + box3.coarseInt();

  ValueType diff = fabs (coarse_int - value);

  if (diff < tol * 1024. ){ // stop
//     value = (value - coarse_int/64.) / (1 - 1./64.);
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




typedef double (*F2) (double, double);

template class Integral2D<F2,double >;
template class Integral2D<F2,float >;




// double indx = 50;
// std::vector<double > global_R246;

// double f (double x, double y)
// {
//   return sin(x) * sin(y);
// }

// double value11_1 (double r, double theta)
// {
//   double COS = cos(theta);
//   double SIN = sin(theta);
//   double COS2 = COS * COS;
//   double COS4 = COS2 * COS2;
//   double SIN2 = SIN * SIN;
//   double SIN4 = SIN2 * SIN2;
//   double r2 = r * r;
//   double r4 = r2 * r2;
//   double R246Q2 = r2 * (global_R246[0] * COS2 + global_R246[1] * COS * SIN + global_R246[2] * SIN2);
//   double R246Q4 = r2 * r2 * (global_R246[0+3] * COS4 + global_R246[1+3] * COS2 * COS * SIN + global_R246[2+3] * COS2 * SIN2 + global_R246[3+3] * COS * SIN * SIN2 + global_R246[4+3] * SIN4);
//   double R246Q6 = r2 * r2 * r2 * (global_R246[0+8] * COS2 * COS2 * COS2 + global_R246[1+8] * COS2 * COS2 * COS * SIN + global_R246[2+8] * COS2 * COS2 * SIN2 + global_R246[3+8] * COS2 * COS * SIN * SIN2 + global_R246[4+8] * COS2 * SIN2 * SIN2 + global_R246[5+8] * COS * SIN * SIN2 * SIN2 + global_R246[6+8] * SIN2 * SIN2 * SIN2);
//   return r4 * r4 * r4 * COS4 * COS4 * COS2 * COS * SIN * r * pow(1-r2,indx) * exp (R246Q2 + R246Q4 + R246Q6);
// }


// double value0_12 (double r, double theta)
// {
//   double COS = cos(theta);
//   double SIN = sin(theta);
//   double COS2 = COS * COS;
//   double COS4 = COS2 * COS2;
//   double SIN2 = SIN * SIN;
//   double SIN4 = SIN2 * SIN2;
//   double r2 = r * r;
//   double r4 = r2 * r2;
//   double R246Q2 = r2 * (global_R246[0] * COS2 + global_R246[1] * COS * SIN + global_R246[2] * SIN2);
//   double R246Q4 = r2 * r2 * (global_R246[0+3] * COS4 + global_R246[1+3] * COS2 * COS * SIN + global_R246[2+3] * COS2 * SIN2 + global_R246[3+3] * COS * SIN * SIN2 + global_R246[4+3] * SIN4);
//   double R246Q6 = r2 * r2 * r2 * (global_R246[0+8] * COS2 * COS2 * COS2 + global_R246[1+8] * COS2 * COS2 * COS * SIN + global_R246[2+8] * COS2 * COS2 * SIN2 + global_R246[3+8] * COS2 * COS * SIN * SIN2 + global_R246[4+8] * COS2 * SIN2 * SIN2 + global_R246[5+8] * COS * SIN * SIN2 * SIN2 + global_R246[6+8] * SIN2 * SIN2 * SIN2);
//   return r4 * r4 * r4 * SIN4 * SIN4 * SIN4 * r * pow(1-r2,indx) * exp (R246Q2 + R246Q4 + R246Q6);
// }

// double value20 (double r, double theta)
// {
//   double COS = cos(theta);
//   double SIN = sin(theta);
//   double COS2 = COS * COS;
//   double SIN2 = SIN * SIN;
//   double r2 = r * r;
//   double R246Q2 = r2 * (global_R246[0] * COS2 + global_R246[1] * COS * SIN + global_R246[2] * SIN2);
//   double R246Q4 = r2 * r2 * (global_R246[0+3] * COS2 * COS2 + global_R246[1+3] * COS2 * COS * SIN + global_R246[2+3] * COS2 * SIN2 + global_R246[3+3] * COS * SIN * SIN2 + global_R246[4+3] * SIN2 * SIN2);
//   double R246Q6 = r2 * r2 * r2 * (global_R246[0+8] * COS2 * COS2 * COS2 + global_R246[1+8] * COS2 * COS2 * COS * SIN + global_R246[2+8] * COS2 * COS2 * SIN2 + global_R246[3+8] * COS2 * COS * SIN * SIN2 + global_R246[4+8] * COS2 * SIN2 * SIN2 + global_R246[5+8] * COS * SIN * SIN2 * SIN2 + global_R246[6+8] * SIN2 * SIN2 * SIN2);
//   return r2 * COS2 * r * pow(1-r2,indx) * exp (R246Q2 + R246Q4 + R246Q6);
// }


// int main ()
// {
//   Integral2D<F2,double > integ;

//   std::vector<double > R246 (15, 0.000);
//   global_R246 = R246;
  
//   double value0;
//   double value1;
//   double e = 2.9615315361867009e-04;
//   // 1.914677683422673e-10;
//   e = 4;
//   double tol = 1e-16;
  
// //   for (int i = 0; i < 1; i ++){
// // //     value0 = integ.cal_int (Integral2DInfo::Gauss9, f, 0, M_PI, 0, M_PI, 1e-10);
// // //     value1 = integ.cal_int (Integral2DInfo::Simpson, f, 0, M_PI, 0, M_PI, 1e-8);
// //     value0 = integ.cal_int (Integral2DInfo::Gauss9, value20, 0, 1, 0, M_PI, tol, 3, 4);
// //     value1 = integ.cal_int (Integral2DInfo::Simpson, value20, 0, 1, 0, M_PI, tol, 3, 4);
// //   }
  
  
//   std::cout.precision (16);
// //   std::cout << std::scientific << value0 - e << std::endl;
// //   std::cout << std::scientific << value1 - e << std::endl;
// //   std::cout << integ.cal_int (Integral2DInfo::Gauss9, value20, 0, 1, 0, M_PI, 1e-17, 3, 4) << std::endl;
  
//   value0 = integ.cal_int (Integral2DInfo::Gauss16, value20, 0, 1, 0, M_PI, tol, 3, 4);
//   value1 = integ.cal_int (Integral2DInfo::Gauss9, value20, 0, 1, 0, M_PI, tol, 3, 4);

//   std::cout << std::scientific << value0  << std::endl;
//   std::cout << std::scientific << value1  << std::endl;
  
//   return 0;
// }
