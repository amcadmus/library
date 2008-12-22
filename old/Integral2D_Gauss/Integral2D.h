
#ifndef __Integral2D_h__wanghan__
#define __Integral2D_h__wanghan__


#include <iostream>
#include <vector>

#define INIT_N 1

typedef double (*F2) (double, double);

template <class BinaryFunction, class ValueType >
class IntegralBox2D;

struct Integral2DMethod 
{
public:
  typedef int Type;
  static Type Simpson;
  static Type Gauss9;
}
    ;

template <class BinaryFunction, class ValueType >
class Integral2D 
{
  ValueType hx, hy;
  IntegralBox2D<BinaryFunction,ValueType >  box;
public:

  ValueType cal_int(const Integral2DMethod::Type method,
		    const BinaryFunction & f, 
		    const ValueType & xlo, const ValueType & xup,
		    const ValueType & ylo, const ValueType & yup,
		    const ValueType & tol = 1e-6, 
		    const int & init_nx = INIT_N,
		    const int & init_ny = INIT_N);
}
    ;



#endif
