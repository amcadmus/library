#ifndef __BisectEqnSolver_h_wanghan__
#define __BisectEqnSolver_h_wanghan__


template <typename VARIABLETYPE, typename VALUETYPE, typename FUNCTION>
class BisectEqnSolver
{
public:
  bool solve (FUNCTION f,
	      const VARIABLETYPE & a,
	      const VARIABLETYPE & b,
	      const VARIABLETYPE & tol,
	      VARIABLETYPE & x);
};





template <typename VARIABLETYPE, typename VALUETYPE, typename FUNCTION>
bool BisectEqnSolver<VARIABLETYPE, VALUETYPE, FUNCTION>::
solve (FUNCTION f,
       const VARIABLETYPE & aa,
       const VARIABLETYPE & bb,
       const VARIABLETYPE & tol,
       VARIABLETYPE & x)
{
  VARIABLETYPE a(aa), b(bb);
  VALUETYPE fa, fb;
  fa = f(a);
  fb = f(b);
  if (fa * fb >0) return false;
  if (fa == VALUETYPE(0)){
    x = a;
    return true;
  }
  if (fb == VALUETYPE(0)){
    x = b;
    return true;
  }
  
  VARIABLETYPE mid;
  VALUETYPE fmid;
  while (fabs(a-b) > tol){
    mid = 0.5 * (a+b);
    fmid = f(mid);
    if (fmid * fa >= 0){
      fa = fmid;
      a = mid;
    }
    else {
      fb = fmid;
      b = mid;
    }
  }

  x = mid;
  return true;
}



#endif
