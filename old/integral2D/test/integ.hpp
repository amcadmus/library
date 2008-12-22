#ifndef _integ_hpp_
#define _integ_hpp_

typedef double (* F_2_DOUBLE) (double, double) ;

double int_1 (F_2_DOUBLE f, const double & a, const double & b, const double & c, const double & d, const double & tol = 1e-6, const int & nx = 20);
double int_1_4 (F_2_DOUBLE f, const double & a, const double & b, const double & c, const double & d, const double & tol = 1e-6, const int & nx = 20);

#endif
