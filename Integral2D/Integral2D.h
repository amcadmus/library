/**
 * @file   Integral2D.h
 * @author Han Wang
 * @date   Mon Jul 30 23:57:59 2007
 * 
 * @brief this is a program to calculate 2-dimensional integration by
 * a self-adaptive integration method. 
 */


#ifndef __Integral2D_h__wanghan__
#define __Integral2D_h__wanghan__

#include <iostream>
#include <vector>

#define INIT_N 1

template <class BinaryFunction, class ValueType >
class IntegralBox2DSimpson;
template <class BinaryFunction, class ValueType >
class IntegralBox2DGauss9;
template <class BinaryFunction, class ValueType >
class IntegralBox2DGauss16;

struct Integral2DInfo 
{
public:
  typedef int Method;
  static Method Simpson;
  static Method Gauss9;
  static Method Gauss16;
}
    ;

template <class BinaryFunction, class ValueType >
class Integral2D 
{
//   BinaryFunction f;
//   ValueType xlo, xup, ylo, yup;
//   ValueType tol;
//   int nx;
//   int ny;
  ValueType hx, hy;
  IntegralBox2DSimpson<BinaryFunction,ValueType >  box_Simpson;
  IntegralBox2DGauss9 <BinaryFunction,ValueType >  box_Gauss9;
  IntegralBox2DGauss16<BinaryFunction,ValueType >  box_Gauss16;
public:

  /** 
   * main function that calculate the integration of a 2-dimensional function.
   *
   * @param method	the integration method
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

  ValueType cal_int (const Integral2DInfo::Method & method,
		     const BinaryFunction & f, 
		     const ValueType & xlo, const ValueType & xup,
		     const ValueType & ylo, const ValueType & yup,
		     const ValueType & tol = 1e-6, 
		     const int & init_nx = INIT_N,
		     const int & init_ny = INIT_N);
}
    ;


template <class BinaryFunction, class ValueType >
class IntegralBox2DSimpson 
{
public:
  IntegralBox2DSimpson () {};

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


template <class BinaryFunction, class ValueType >
class IntegralBox2DGauss9 
{
public:
  IntegralBox2DGauss9 () {};

  void reinit (const BinaryFunction & f, 
	       const ValueType & xlo, const ValueType & xup,
	       const ValueType & ylo, const ValueType & yup,
	       const ValueType & tol );
  ValueType & coarseInt () {return coarse_int;}
  const ValueType & coarseInt () const {return coarse_int;}
  
  void cal_int (ValueType & value);
private:
  BinaryFunction f;
  ValueType hx,  hy;
  ValueType xlo, xup;
  ValueType ylo, yup;
  std::vector<ValueType > coarse_value;
  ValueType coarse_int;
  ValueType tol;
}
    ;

template <class BinaryFunction, class ValueType >
class IntegralBox2DGauss16 
{
public:
  IntegralBox2DGauss16 () {};

  void reinit (const BinaryFunction & f, 
	       const ValueType & xlo, const ValueType & xup,
	       const ValueType & ylo, const ValueType & yup,
	       const ValueType & tol );
  ValueType & coarseInt () {return coarse_int;}
  const ValueType & coarseInt () const {return coarse_int;}
  
  void cal_int (ValueType & value);
private:
  BinaryFunction f;
  ValueType hx,  hy;
  ValueType xlo, xup;
  ValueType ylo, yup;
  std::vector<ValueType > coarse_value;
  ValueType coarse_int;
  ValueType tol;
}
    ;


#endif	// __Integral2D_h__wanghan__

/**
 * end of file
 * 
 */
