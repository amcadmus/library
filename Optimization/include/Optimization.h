#ifndef __Optimization_h_wanghan__
#define __Optimization_h_wanghan__

#include <vector>
#include <stdio.h>
#include "Target.h"

#define  ScalorType double

// class OPTIMIZATION_TARGET 
// {
//   unsigned natom;
//   const std::vector<std::vector<double > > &posis;
//   double mode2;
//   double mylambda;
// public:
//   OPTIMIZATION_TARGET (const std::vector<std::vector<double > > & p_posis,
// 		       const ScalorType & lambda);
//   void value (const std::vector<ScalorType > & x,
// 	      ScalorType & value) const;
//   void grad  (const std::vector<ScalorType > & x,
// 	      std::vector<ScalorType > & g) const;
//   void hessian (const std::vector<ScalorType > & x,
// 		std::vector<ScalorType > & hessian) const;
// };

void
optimization_printVec (const std::vector<ScalorType> & x);


// class OPTIMIZATION_TARGET
// {
//   double hg, hh;
//   unsigned natom;
//   const std::vector<std::vector<double > > &posis;
//   double mode2;
//   double mylambda;
// public:
//   OPTIMIZATION_TARGET (const std::vector<std::vector<double > > & p_posis,
// 			const ScalorType & lambda);
//   void value (const std::vector<ScalorType > & x,
// 	      ScalorType & value) const;
//   void grad  (const std::vector<ScalorType > & x,
// 	      std::vector<ScalorType > & g) const;
//   void hessian (const std::vector<ScalorType > & x,
// 		std::vector<ScalorType > & hessian) const;
// };



class OptimizationCounter 
{
  unsigned cv;
  unsigned cg;
  unsigned ch;
  unsigned nstep;
public:
  OptimizationCounter ();
  void reset ();
  const unsigned & count_value   (const unsigned & c);
  const unsigned & count_grad    (const unsigned & c);
  const unsigned & count_hessian (const unsigned & c);
  const unsigned & count_step    ();
  const unsigned & getCount_value   () const {return cv;}
  const unsigned & getCount_grad    () const {return cg;}
  const unsigned & getCount_hessian () const {return ch;}
  const unsigned & getCount_step    () const {return nstep;}
  void print () const;
};


struct OptimizationControler 
{
  unsigned maxIter;
  bool printDetail;
  ScalorType desiredPrecision_stepDiff;
  ScalorType desiredPrecision_gradNorm;
  OptimizationControler ();
};


class LineSearch 
{
public:
  ScalorType exactLineSearchDefaultPrecision;
  unsigned searchCubicMaxIter;
  unsigned searchBisectMaxIter;
  ScalorType WP_sigma;
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
  bool exactSearch_618 (const OPTIMIZATION_TARGET & f,
			const std::vector<ScalorType > & x,
			const std::vector<ScalorType > & d,
			const ScalorType & delta,
			ScalorType & c,
			OptimizationCounter * oc = NULL);
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



#endif
