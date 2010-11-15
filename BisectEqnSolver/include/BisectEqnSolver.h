#ifndef __BisectEqnSolver_h_wanghan__
#define __BisectEqnSolver_h_wanghan__

typedef  double VARIABLETYPE;
typedef double VALUETYPE;
typedef VALUETYPE (*PDF) (VALUETYPE);

class BisectEqnSolver
{
public:
  bool solve (PDF f,
	      const VARIABLETYPE & a,
	      const VARIABLETYPE & b,
	      const VARIABLETYPE & tol,
	      VARIABLETYPE & x);
};


#endif
