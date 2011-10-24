#include "Optimization.h"
#include "CholeskySolver.h"

#include <cmath>

void
optimization_printVec (const std::vector<ScalorType> & x)
{
  for (unsigned i = 0; i < x.size(); ++i){
    if (x[i]>=0) putchar(' ');
    printf ("%.8e  ", x[i]);
  }
  for (unsigned i = 0; i < 3; ++i) putchar (' ');
}

OptimizationControler::
OptimizationControler ()
    : maxIter (-1),
      printDetail (false),
      desiredPrecision_stepDiff(1e-4),
      desiredPrecision_gradNorm(1e-4)
{
}

OptimizationCounter::
OptimizationCounter ()
    : cv(0), cg(0), ch(0), nstep(0)
{
}

inline void OptimizationCounter::
reset ()
{
  cv = cg = ch = 0;
}

inline const unsigned & OptimizationCounter::
count_value (const unsigned & c)
{
  return cv += c;
}

inline const unsigned & OptimizationCounter::
count_grad (const unsigned & c)
{
  return cg += c;
}

inline const unsigned & OptimizationCounter::
count_hessian (const unsigned & c)
{
  return ch += c;
}

inline const unsigned & OptimizationCounter::
count_step ()
{
  return ++ nstep;
}

void OptimizationCounter::
print () const
{
  printf ("# number of iter    is %d\n", nstep);
  printf ("# number of value   is %d\n", cv);
  printf ("# number of frad    is %d\n", cg);
  printf ("# number of hessian is %d\n", ch);
}


LineSearch::
LineSearch()
    : exactLineSearchDefaultPrecision(1e-2),
      searchCubicMaxIter(30),
      searchBisectMaxIter(30),
      WP_sigma(0.1)
{
}


void LineSearch::
initializeZone (const OPTIMIZATION_TARGET & f,
		const std::vector<ScalorType > & x,
		const std::vector<ScalorType > & d,
		ScalorType & a,
		ScalorType & b,
		OptimizationCounter * oc)
{
  a = 0;
  ScalorType h = 1, t = 3;
  b = a + h;
  ScalorType fa, fb;
  std::vector<ScalorType > ad(x), bd(x);
  for (unsigned i = 0; i < x.size(); ++i){
    ad[i] += a * d[i];
    bd[i] += b * d[i];
  }
  f.value (ad, fa);
  f.value (bd, fb);
  if (oc != NULL) oc->count_value(2);
  
  while (fa > fb){
    h *= t;
    a = b;
    fa = fb;
    b = a + h;
    for (unsigned i = 0; i < x.size(); ++i){
      bd[i] = x[i] + b * d[i];
    }
    f.value (bd, fb);
    if (oc != NULL) oc->count_value(1);
  }
  
  a = 0;
}

ScalorType LineSearch::
phiP (const OPTIMIZATION_TARGET & f,
      const std::vector<ScalorType > & x,
      const std::vector<ScalorType > & d,
      const ScalorType & a,
      OptimizationCounter * oc)
{
  std::vector<ScalorType > x1 (x);
  for (unsigned i = 0; i < x.size(); ++i){
    x1[i] += a * d[i];
  }
  std::vector<ScalorType > g;
  f.grad (x1, g);
  if (oc != NULL) oc->count_grad(1);
  ScalorType sum = ScalorType(0.);
  for (unsigned i = 0; i < x.size(); ++i){
    sum += g[i] * d[i];
  }
  return sum;
}

ScalorType LineSearch::
phi (const OPTIMIZATION_TARGET & f,
     const std::vector<ScalorType > & x,
     const std::vector<ScalorType > & d,
     const ScalorType & a,
     OptimizationCounter * oc)
{
  std::vector<ScalorType > x1 (x);
  for (unsigned i = 0; i < x.size(); ++i){
    x1[i] += a * d[i];
  }
  ScalorType fa = 0.;
  f.value (x1, fa);
  if (oc != NULL) oc->count_value(1);
  
  return fa;
}


bool LineSearch::
softSearch_cubic (const OPTIMIZATION_TARGET & f,
		  const std::vector<ScalorType > & x,
		  const std::vector<ScalorType > & d,
		  ScalorType & c,
		  OptimizationCounter * oc)
{
  ScalorType a, b;
  initializeZone (f, x, d, a, b, oc);

  ScalorType phia, phib, phipa, phipb, phip0;
  phia = phi (f, x, d, a, oc);
  phib = phi (f, x, d, b, oc);
  phipa = phiP (f, x, d, a, oc);
  phipb = phiP (f, x, d, b, oc);
  phip0 = phiP (f, x, d, ScalorType(0.), oc);

  unsigned count = 0;
  while (true){  
    ScalorType s = 3. * (phib - phia) / (b - a);
    ScalorType z = s - phipb - phipa;
    ScalorType w = z * z - phipa * phipb;
    if (w >= 0){
      w = sqrt(w);
      c = a + (b - a) * (1 - (w - phipa - z) / (phipb - phipa + 2. * w));
    }
    else {
      c = 0.5 * (a + b);
    }

    ScalorType phipc = phiP (f, x, d, c, oc);
    if (phipc < - WP_sigma * phip0){
      if (phipc  > WP_sigma * phip0){
	return true;
      }
      else {
	if ((count ++) == searchCubicMaxIter){
	  break;
	}
	a = c;
	phia = phi(f, x, d, a, oc);
	phipa = phipc;
      }
    }
    else {
      if ((count ++) == searchCubicMaxIter){
	break;
      }
      b = c;
      phib = phi(f, x, d, b, oc);
      phipb = phipc;
    }
  }
  return false;
}


bool LineSearch::
softSearch_bisect (const OPTIMIZATION_TARGET & f,
		   const std::vector<ScalorType > & x,
		   const std::vector<ScalorType > & d,
		   ScalorType & c,
		   OptimizationCounter * oc)
{
  ScalorType a, b;
  initializeZone (f, x, d, a, b, oc);

  ScalorType phip0;
  phip0 = phiP (f, x, d, ScalorType(0.), oc);

  unsigned count = 0;
  while (true){  
    c = 0.5 * (a + b);
    ScalorType phipc = phiP (f, x, d, c, oc);
    
    if (phipc < - WP_sigma * phip0){
      if (phipc  > WP_sigma * phip0){
	return true;
      }
      else {
	if ((count ++) == searchBisectMaxIter){
	  break;
	}
	a = c;
      }
    }
    else {
      if ((count ++) == searchBisectMaxIter){
	break;
      }
      b = c;
    }
  }
  return false;
}


bool LineSearch::
exactSearch_618 (const OPTIMIZATION_TARGET & f,
		 const std::vector<ScalorType > & x,
		 const std::vector<ScalorType > & d,
		 const ScalorType & delta,
		 ScalorType & c,
		 OptimizationCounter * oc)
{  
  ScalorType a, b;
  initializeZone (f, x, d, a, b, oc);
  ScalorType lambda = a + 0.382 * (b - a);
  ScalorType mu     = a + 0.618 * (b - a);
  ScalorType phil = phi(f, x, d, lambda, oc);
  ScalorType phim = phi(f, x, d, mu,     oc);

  while (true){
    if (phil > phim){
      if (b - lambda <= delta) {
	c = mu;
	return true;
      }
      else{
	a = lambda;
	lambda = mu;
	phil = phim;
	mu = a + 0.618 * (b - a);
	phim = phi(f, x, d, mu,     oc);
      }
    }
    else {
      if (mu - a <= delta){
	c = lambda;
	return true;
      }
      else {
	b = mu;
	mu = lambda;
	phim = phil;
	lambda = a + 0.382 * (b - a);
	phil = phi(f, x, d, lambda, oc);
      }
    }
  }
  
  return false;
}



bool LineSearch::
softSearch_mix (const OPTIMIZATION_TARGET & f,
		const std::vector<ScalorType > & x,
		const std::vector<ScalorType > & d,
		ScalorType & c,
		OptimizationCounter * oc)
{
  if (softSearch_cubic(f, x, d, c, oc)){
    return true;
  }
  else {
    if (softSearch_bisect (f, x, d, c, oc)){
      return true;
    }
    else{
      return false;
    }
  }
}  


void SpeedestDescent::
findMin (const OPTIMIZATION_TARGET & f,
	 const OptimizationControler & control,
	 std::vector<ScalorType > & x,
	 OptimizationCounter * oc)
{
  LineSearch ls;
  unsigned step = 0;
  
  std::vector<ScalorType > g;
  f.grad (x, g);
  if (oc != NULL) oc->count_grad(1);
  ScalorType norm = 0;
  for (unsigned i = 0; i < x.size(); ++i){
    g[i] *= -1.;
    norm += g[i] * g[i];
  }
  ScalorType e2 = control.desiredPrecision_gradNorm;
  e2 *= e2;

  if (control.printDetail) {
    printf ("# step %08d    x = ", step);
    optimization_printVec (x);
    ScalorType tmpvalue;
    f.value(x, tmpvalue);
    if (oc != NULL) oc->count_value(1);
    printf ("f = %e, ||g||^2 = %e, e^2 = %e\n", tmpvalue, norm, e2);
  }
  
  while (norm > e2 && step < control.maxIter){
    ScalorType c;
    if (! ls.softSearch_mix (f, x, g, c, oc)){
      ScalorType delta = control.desiredPrecision_stepDiff  / norm;
      if (delta > ls.exactLineSearchDefaultPrecision) {
	delta = ls.exactLineSearchDefaultPrecision;
      }
      ls.exactSearch_618 (f, x, g, delta, c, oc);
    }
    for (unsigned i = 0; i < x.size(); ++i){
      x[i] += c * g[i];
    }
    f.grad (x, g);
    if (oc != NULL) oc->count_grad(1);
    norm = 0.;
    for (unsigned i = 0; i < x.size(); ++i){
      g[i] *= -1.;
      norm += g[i] * g[i];
    }
    step ++;
    if (control.printDetail) {
      printf ("# step %08d    x = ", step);
      optimization_printVec (x);
      ScalorType tmpvalue;
      f.value(x, tmpvalue);
      if (oc != NULL) oc->count_value(1);
      printf ("f = %e, ||g||^2 = %e, e^2 = %e\n", tmpvalue, norm, e2);
    }    
    if (oc != NULL) oc->count_step();
  }
}


void NewtonInteration::
findMin (const OPTIMIZATION_TARGET & f,
	 const OptimizationControler & control,
	 std::vector<ScalorType > & x,
	 OptimizationCounter * oc)
{
  CholeskySolver<ScalorType> cs;
  LineSearch ls;  
  std::vector<ScalorType > g, H;
  unsigned step = 0;
  
  f.grad (x, g);
  if (oc != NULL) oc->count_grad(1);
  ScalorType norm = 0;
  for (unsigned i = 0; i < x.size(); ++i){
    g[i] *= -1.;
    norm += g[i] * g[i];
  }
  ScalorType e2 = control.desiredPrecision_gradNorm;
  e2 *= e2;
  
  if (control.printDetail) {
    printf ("# step %08d    x = ", step);
    optimization_printVec (x);
    ScalorType tmpvalue;
    f.value(x, tmpvalue);
    if (oc != NULL) oc->count_value(1);
    for (unsigned i = 0; i < 13; ++i) putchar(' ');
    printf ("f = %e, ||g||^2 = %e, e^2 = %e\n", tmpvalue, norm, e2);
  }
  unsigned type = 0; // 0: Newton step, 1: SD step
  
  while (norm > e2 && step < control.maxIter){
    f.hessian (x, H);
    if (oc != NULL) oc->count_hessian(1);
    if (cs.reinit (H, x.size()) == 0 && cs.solve(g) == 0){
      ScalorType c;
      if (! ls.softSearch_mix (f, x, g, c, oc)){
	ScalorType delta = control.desiredPrecision_stepDiff  / norm;
	if (delta > ls.exactLineSearchDefaultPrecision) {
	  delta = ls.exactLineSearchDefaultPrecision;
	}
	ls.exactSearch_618 (f, x, g, delta, c, oc);
      }
      for (unsigned i = 0; i < x.size(); ++i){
	x[i] += c * g[i];
      }
      type = 0;
    }
    else {
      ScalorType c;
      ls.softSearch_mix (f, x, g, c, oc);
      for (unsigned i = 0; i < x.size(); ++i){
	x[i] += c * g[i];
      }
      type = 1;
    }
    
    f.grad (x, g);
    if (oc != NULL) oc->count_grad(1);
    norm = 0.;
    for (unsigned i = 0; i < x.size(); ++i){
      g[i] *= -1.;
      norm += g[i] * g[i];
    }
    step ++;
    
    if (control.printDetail) {
      printf ("# step %08d    x = ", step);
      optimization_printVec (x);
      ScalorType tmpvalue;
      f.value(x, tmpvalue);
      if (oc != NULL) oc->count_value(1);
      if (type == 0) printf ("Newton step: ");
      if (type == 1) printf ("SD     step: ");      
      printf ("f = %e, ||g||^2 = %e, e^2 = %e\n", tmpvalue, norm, e2);
    }
    if (oc != NULL) oc->count_step();
  }
}

    
