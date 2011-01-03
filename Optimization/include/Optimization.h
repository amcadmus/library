#ifndef __Optimization_h_wanghan__
#define __Optimization_h_wanghan__

#include <vector>
#include <stdio.h>

#include "Target.h"

typedef  double ScalorType;

class OptimizationCounter 
{
  unsigned cv;
  unsigned cg;
  unsigned ch;
public:
  OptimizationCounter ();
  void reset ();
  const unsigned & count_value   (const unsigned & c);
  const unsigned & count_grad    (const unsigned & c);
  const unsigned & count_hessian (const unsigned & c);
  const unsigned & getCount_value   () const {return cv;}
  const unsigned & getCount_grad    () const {return cg;}
  const unsigned & getCount_hessian () const {return ch;}
  void print () const;
};


struct OptimizationControler 
{
  unsigned maxIter;
  ScalorType desiredPrecision_stepDiff;
  ScalorType desiredPrecision_gradNorm;
  OptimizationControler ();
};



class LineSearch 
{
  unsigned searchCubicMaxIter;
  unsigned searchBisectMaxIter;
  ScalorType sigma;
private:
  void initializeZone (const OPTIMIZATION_TARGET & f,
		       const std::vector<ScalorType > & x,
		       const std::vector<ScalorType > & d,
		       ScalorType & a,
		       ScalorType & b,
		       OptimizationCounter * oc = NULL);
  ScalorType phi  (const OPTIMIZATION_TARGET & f,
		   const std::vector<ScalorType > & x,
		   const std::vector<ScalorType > & d,
		   const ScalorType & a,
		   OptimizationCounter * oc = NULL);
  ScalorType phiP (const OPTIMIZATION_TARGET & f,
		   const std::vector<ScalorType > & x,
		   const std::vector<ScalorType > & d,
		   const ScalorType & a,
		   OptimizationCounter * oc = NULL);
public:
  LineSearch ();
  bool softSearch_cubic (const OPTIMIZATION_TARGET & f,
			 const std::vector<ScalorType > & x,
			 const std::vector<ScalorType > & d,
			 ScalorType & c,
			 OptimizationCounter * oc = NULL);
  bool softSearch_bisect (const OPTIMIZATION_TARGET & f,
			  const std::vector<ScalorType > & x,
			  const std::vector<ScalorType > & d,
			  ScalorType & c,
			  OptimizationCounter * oc = NULL);
  bool softSearch_mix (const OPTIMIZATION_TARGET & f,
		       const std::vector<ScalorType > & x,
		       const std::vector<ScalorType > & d,
		       ScalorType & c,
		       OptimizationCounter * oc = NULL);
  
};



class SpeedestDescent 
{
public:
  void findMin (const OPTIMIZATION_TARGET & f,
		const OptimizationControler & control,
		std::vector<ScalorType > & x,
		OptimizationCounter * oc);
};  


class NewtonInteration
{
public:
  void findMin (const OPTIMIZATION_TARGET & f,
		const OptimizationControler & control,
		std::vector<ScalorType > & x,
		OptimizationCounter * oc);
};


static void
printVec (const std::vector<ScalorType> & x)
{
  for (unsigned i = 0; i < x.size(); ++i){
    printf ("%f  ", x[i]);
  }
  printf ("\n");
}

#endif
