/**
 * @file   Integral2D.h
 * @author Han Wang
 * @date   Mon Jul 30 23:57:59 2007
 * 
 * @brief this is a program to calculate 2-dimensional integration by
 * a self-adaptive Simpson method. 
 */


#ifndef __Integral2D_h__wanghan__
#define __Integral2D_h__wanghan__

#include <iostream>
#include <vector>

#define INIT_N 1

template <class BinaryFunction, class ValueType >
class IntegralBox2D;

template <class BinaryFunction, class ValueType >
class Integral2D 
{
//   BinaryFunction f;
//   ValueType xlo, xup, ylo, yup;
//   ValueType tol;
//   int nx;
//   int ny;
  ValueType hx, hy;
  IntegralBox2D<BinaryFunction,ValueType >  box;
public:
  /** 
   * main function that calculate the integration of a 2-dimensional function.
   * 
   * @param f		the name of a binary function being integrated. it will be call as ValueType value = f (ValueType, ValueType)
   * @param xlo		the lower boundary of the first dimension
   * @param xup		the upper boundary of the first dimension
   * @param ylo		the lower boundary of the second dimension
   * @param yup		the upper boundary of the second dimension
   * @param tol		absolute tolerance
   * @param init_nx	number of cells on the first dimension of the initial mesh. default value is 1
   * @param init_ny	number of cells on the second dimension of the initial mesh. default value is 1
   * 
   * @return		the integration value
   */

  ValueType cal_int(const BinaryFunction & f, 
		    const ValueType & xlo, const ValueType & xup,
		    const ValueType & ylo, const ValueType & yup,
		    const ValueType & tol = 1e-6, 
		    const int & init_nx = INIT_N,
		    const int & init_ny = INIT_N);
}
    ;


template <class BinaryFunction, class ValueType >
class IntegralBox2D 
{
public:
  IntegralBox2D () {coarse_value.resize (9);};

/** 
 * initialize function.
 * 
 * @param f		the binary function to be integrated. it will be call as ValueType value = f (ValueType, ValueType)
 * @param xlo		the lower boundary of the first dimension of the box
 * @param xup		the upper boundary of the first dimension of the box
 * @param ylo		the lower boundary of the second dimension of the box
 * @param yup		the upper boundary of the second dimension of the box
 * @param tol		absolute tolerance
 */
  void reinit (const BinaryFunction & f, 
	       const ValueType & xlo, const ValueType & xup,
	       const ValueType & ylo, const ValueType & yup,
	       const ValueType & tol );
/** 
 * initialize function. if the result of 9 points formula (Simpson
 * integration formula) on the box is known
 * 
 * @param f		the binary function to be integrated. it will be call as ValueType value = f (ValueType, ValueType)
 * @param xlo		the lower boundary of the first dimension of the box
 * @param xup		the upper boundary of the first dimension of the box
 * @param ylo		the lower boundary of the second dimension of the box
 * @param yup		the upper boundary of the second dimension of the box
 * @param coarse_value	result of 9 points formula on the box
 * @param coarse_int	value of the 9 integration points
 * @param tol		absolute tolerance
 */
  void reinit (const BinaryFunction & f, 
	       const ValueType & xlo, const ValueType & xup, 
	       const ValueType & ylo, const ValueType & yup,
	       const std::vector<ValueType > & coarse_value,
	       const ValueType & coarse_int,
	       const ValueType & tol);
/** 
 * calculate the integration of f. first the box is divided into four
 * refined boxes, then 9 points formula is acted on them, at last the
 * error is estimated by the refined values and coarse value to decide
 * whether the box should be refined.
 * 
 * @param value		the integration of f
 */
  void cal_int (ValueType & value);
private:
  BinaryFunction f;
  ValueType hx,  hy;
  ValueType xlo, xup;
  ValueType ylo, yup;
  ValueType one_36;
  std::vector<ValueType > coarse_value;
  ValueType coarse_int;
  ValueType tol;
}
    ;

////////////////////////////////////////////////////////////////////////////////////////////////////

static double global_hx;
static double global_hy;

template <class BinaryFunction, class ValueType >
ValueType Integral2D<BinaryFunction,ValueType >::cal_int (const BinaryFunction & f, 
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
  hx = 0.5 * (xup - xlo);
  hy = 0.5 * (yup - ylo);
  one_36 = 1./36.;
  
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
void IntegralBox2D<BinaryFunction,ValueType >::reinit (const BinaryFunction & f_, 
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
void IntegralBox2D<BinaryFunction,ValueType >::cal_int (ValueType & value)
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

  if (diff < tol * 15. ){ // stop
    return;
  }
  else{
    if ( hx < 1e-14 * global_hx || hy < 1e-14 * global_hy){
//     std::cout << "cannot reach the tolrence. " << "hx " << hx << " hy " << hy << std::endl;
      return;
    }
    else{
      IntegralBox2D box0;
      IntegralBox2D box1;
      IntegralBox2D box2;
      IntegralBox2D box3;
    
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
}

#endif	// __Integral2D_h__wanghan__

/**
 * end of file
 * 
 */
