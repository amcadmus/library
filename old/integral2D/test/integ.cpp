#include<iostream>
#include<math.h>
#include<vector>
#include<list>
#include"integ.hpp"

using namespace std;
double int_12 (F_2_DOUBLE f, const double & rho, const double & c, const double & d, const double & tol = 1e-6, const int & ny = 20);
double int_12_4 (F_2_DOUBLE f, const double & rho, const double & c, const double & d, const double & tol = 1e-6, const int & ny = 20);


double int_1 (F_2_DOUBLE f, const double & a, const double & b, const double & c, const double & d, const double & tol, const int & nx)
{
    double integration = 0;
    
    list<double> r (nx+1);
    list<double> fr (nx+1);
    list<int> tag (nx+1, 0);
    double ddx = (b - a)/nx;
    
    list<double>::iterator i;
    int count = 0;
    for (i = r.begin(); i != r.end(); i++){
	*i = (count ++) * ddx + a;
    }

    list<double>::iterator this_r = r.begin();
    list<double>::iterator next_r = ++ r.begin();
    list<double>::iterator this_fr = fr.begin();
    list<double>::iterator next_fr = ++ fr.begin();
    list<int>::iterator this_tag = tag.begin();
    list<int>::iterator next_tag = ++ tag.begin();
    
    int goon = 0;
    while (next_r != r.end()){
	double dx = * next_r - * this_r;
	if (fabs(dx) < 1e-10){
	    fprintf(stderr, "cannot reach the required tolerance!\n");
	    goon = 1;
	}
	double my_tol = tol * dx / fabs(b-a);
	double point_tol = 0.5 * my_tol / dx;
	double mid_r = 0.5 * (* next_r + * this_r);
	double mid_fr;
	
//	point_tol = 1e-0;
	
	if (* this_tag == 0){
	    * this_fr = int_12 (f, * this_r, c, d, point_tol);
	    * this_tag = 1;
	}
	if (* next_tag == 0){
	    * next_fr = int_12 (f, * next_r, c, d, point_tol);
	    * next_tag = 1;
	}

	mid_fr = int_12 (f, mid_r, c, d, point_tol);
	
// 	double s1 = 0.5 * (* this_fr + * next_fr) * dx;
// 	double s21 = 0.5 * (* this_fr + mid_fr) * 0.5 * dx;
// 	double s22 = 0.5 * (* next_fr + mid_fr) * 0.5 * dx;
	double s1 = 1./6. * (*this_fr + 4*mid_fr + *next_fr) * dx;
	double s21 = 1./6. * (*this_fr + 4*int_12 (f, 0.5*(mid_r+*this_r), c, d, point_tol) + mid_fr) * 0.5 * dx;
	double s22 = 1./6. * (*next_fr + 4*int_12 (f, 0.5*(mid_r+*next_r), c, d, point_tol) + mid_fr) * 0.5 * dx;
	
	if (fabs (s1 - s21 - s22) < 31 * my_tol || goon == 1){
	    integration += s21 + s22;
	    this_r ++;
	    next_r ++;
	    this_fr ++;
	    next_fr ++;
	    this_tag ++;
	    next_tag ++;
	    if (goon == 1){
		goon = 0;
	    }
	}
	else {
	    next_r = r.insert (next_r, mid_r);
	    next_fr = fr.insert (next_fr, mid_fr);
	    next_tag = tag.insert (next_tag, 1);
	}
    }
//    printf("tol=%e\tn_point=%d\n",tol,r.size());
    
    return integration;
}


double int_1_4 (F_2_DOUBLE f, const double & a, const double & b, const double & c, const double & d, const double & tol, const int & nx)
{
    double integration = 0;
    
    list<double> r (nx+1);
    list<double> fr (nx+1);
    list<int> tag (nx+1, 0);
    double ddx = (b - a)/nx;
    
    list<double>::iterator i;
    int count = 0;
    for (i = r.begin(); i != r.end(); i++){
	*i = (count ++) * ddx + a;
    }

    list<double>::iterator this_r = r.begin();
    list<double>::iterator next_r = ++ r.begin();
    list<double>::iterator this_fr = fr.begin();
    list<double>::iterator next_fr = ++ fr.begin();
    list<int>::iterator this_tag = tag.begin();
    list<int>::iterator next_tag = ++ tag.begin();
    
    while (next_r != r.end()){
	double dx = * next_r - * this_r;
	double my_tol = tol * dx / fabs(b-a);
	double point_tol = 0.5 * my_tol / dx;
	double mid1_r = * this_r + 0.25 * dx;
	double mid2_r = * this_r + 0.50 * dx;
	double mid3_r = * this_r + 0.75 * dx;
	
//	point_tol = 1e-0;
	
	if (* this_tag == 0){
	    * this_fr = int_12_4 (f, * this_r, c, d, point_tol);
	    * this_tag = 1;
	}
	if (* next_tag == 0){
	    * next_fr = int_12_4 (f, * next_r, c, d, point_tol);
	    * next_tag = 1;
	}

	double mid1_fr = int_12_4 (f, mid1_r, c, d, point_tol);
	double mid2_fr = int_12_4 (f, mid2_r, c, d, point_tol);
	double mid3_fr = int_12_4 (f, mid3_r, c, d, point_tol);
	
	double s1 = 1./6. * (*this_fr + 4 * mid2_fr + *next_fr) * dx;
	double s21 = 1./6. * (*this_fr + 4 * mid1_fr + mid2_fr) * 0.5 * dx;
	double s22 = 1./6. * (*next_fr + 4 * mid3_fr + mid2_fr) * 0.5 * dx;
	
	if (fabs (s1 - s21 - s22) < 31 * my_tol){
	    integration += s21 + s22;
	    this_r ++;
	    next_r ++;
	    this_fr ++;
	    next_fr ++;
	    this_tag ++;
	    next_tag ++;
	}
	else {
	    next_r = r.insert (next_r, mid3_r);
	    next_r = r.insert (next_r, mid2_r);
	    next_r = r.insert (next_r, mid1_r);
	    next_fr = fr.insert (next_fr, mid3_fr);
	    next_fr = fr.insert (next_fr, mid2_fr);
	    next_fr = fr.insert (next_fr, mid1_fr);
	    next_tag = tag.insert (next_tag, 1);
	    next_tag = tag.insert (next_tag, 1);
	    next_tag = tag.insert (next_tag, 1);
	}
    }
    printf("tol=%e\tn_point=%d\n",tol,r.size());
    
    return integration;
}


double int_12_4 (F_2_DOUBLE f, const double & rho, const double & c, const double & d, const double & tol, const int & ny)
{
    double integration = 0;
    
    list<double> r (ny+1);
    list<double> fr (ny+1);
    list<int> tag (ny+1, 0);
    double ddy = (d - c)/ny;
    
    list<double>::iterator i;
    int count = 0;
    for (i = r.begin(); i != r.end(); i++){
	*i = (count ++) * ddy + c;
    }

    list<double>::iterator this_r = r.begin();
    list<double>::iterator next_r = ++ r.begin();
    list<double>::iterator this_fr = fr.begin();
    list<double>::iterator next_fr = ++ fr.begin();
    list<int>::iterator this_tag = tag.begin();
    list<int>::iterator next_tag = ++ tag.begin();
    
    int goon = 0;
    while (next_r != r.end()){
	double dx = * next_r - * this_r;
	if (fabs(dx) < 1e-13){
	    fprintf(stderr, "cannot reach the required tolerance!\n");
	    goon = 1;
	}
	double my_tol = tol * dx / fabs(d-c);
//	double point_tol = 0.5 * my_tol / dx;
	double mid1_r = * this_r + 0.25 * dx;
	double mid2_r = * this_r + 0.50 * dx;
	double mid3_r = * this_r + 0.75 * dx;
	
	if (* this_tag == 0){
	    * this_fr = f(rho, * this_r);
	    * this_tag = 1;
	}
	if (* next_tag == 0){
	    * next_fr = f(rho, * next_r);
	    * next_tag = 1;
	}
//////////////////////////////
// here maybe optimazed !!
//////////////////////////////
	double mid1_fr = f(rho, mid1_r);
	double mid2_fr = f(rho, mid2_r);
	double mid3_fr = f(rho, mid3_r);
	
// 	double s1 = 0.5 * (* this_fr + * next_fr) * dx;
// 	double s21 = 0.5 * (* this_fr + mid_fr) * 0.5 * dx;
// 	double s22 = 0.5 * (* next_fr + mid_fr) * 0.5 * dx;
	
	double s1 = 1./6. * (*this_fr + 4 * mid2_fr + *next_fr) * dx;
	double s21 = 1./6. * (*this_fr + 4 * mid1_fr + mid2_fr) * 0.5 * dx;
	double s22 = 1./6. * (*next_fr + 4 * mid3_fr + mid2_fr) * 0.5 * dx;

	
	if (fabs (s1 - s21 - s22) < 31 * my_tol || goon == 1){
	    integration += s21 + s22;
	    this_r ++;
	    next_r ++;
	    this_fr ++;
	    next_fr ++;
	    this_tag ++;
	    next_tag ++;
	    if (goon == 1){
		goon = 0;
	    }
	}
	else {
	    next_r = r.insert (next_r, mid3_r);
	    next_r = r.insert (next_r, mid2_r);
	    next_r = r.insert (next_r, mid1_r);
	    next_fr = fr.insert (next_fr, mid3_fr);
	    next_fr = fr.insert (next_fr, mid2_fr);
	    next_fr = fr.insert (next_fr, mid1_fr);
	    next_tag = tag.insert (next_tag, 1);
	    next_tag = tag.insert (next_tag, 1);
	    next_tag = tag.insert (next_tag, 1);
	}
    }
//    printf("tol=%e\tn_point=%d\n",tol,r.size());
    
    return integration;
}



double int_12 (F_2_DOUBLE f, const double & rho, const double & c, const double & d, const double & tol, const int & ny)
{
    double integration = 0;
    
    list<double> r (ny+1);
    list<double> fr (ny+1);
    list<int> tag (ny+1, 0);
    double ddy = (d - c)/ny;
    
    list<double>::iterator i;
    int count = 0;
    for (i = r.begin(); i != r.end(); i++){
	*i = (count ++) * ddy + c;
    }

    list<double>::iterator this_r = r.begin();
    list<double>::iterator next_r = ++ r.begin();
    list<double>::iterator this_fr = fr.begin();
    list<double>::iterator next_fr = ++ fr.begin();
    list<int>::iterator this_tag = tag.begin();
    list<int>::iterator next_tag = ++ tag.begin();
    
    int goon = 0;
    while (next_r != r.end()){
	double dx = * next_r - * this_r;
	if (fabs(dx) < 1e-13){
	    fprintf(stderr, "cannot reach the required tolerance!\n");
	    goon = 1;
	}
	double my_tol = tol * dx / fabs(d-c);
//	double point_tol = 0.5 * my_tol / dx;
	double mid_r = 0.5 * (* next_r + * this_r);
	double mid_fr;
	
	if (* this_tag == 0){
	    * this_fr = f(rho, * this_r);
	    * this_tag = 1;
	}
	if (* next_tag == 0){
	    * next_fr = f(rho, * next_r);
	    * next_tag = 1;
	}
//////////////////////////////
// here maybe optimazed !!
//////////////////////////////
	mid_fr = f(rho, mid_r);
	
// 	double s1 = 0.5 * (* this_fr + * next_fr) * dx;
// 	double s21 = 0.5 * (* this_fr + mid_fr) * 0.5 * dx;
// 	double s22 = 0.5 * (* next_fr + mid_fr) * 0.5 * dx;
	double s1 = 1./6. * (*this_fr + 4*mid_fr + *next_fr) * dx;
	double s21 = 1./6. * (*this_fr + 4 * f(rho, 0.5*(mid_r+*this_r)) + mid_fr) * 0.5 * dx;
	double s22 = 1./6. * (*next_fr + 4 * f(rho, 0.5*(mid_r+*next_r)) + mid_fr) * 0.5 * dx;

	
	if (fabs (s1 - s21 - s22) < 31 * my_tol || goon == 1){
	    integration += s21 + s22;
	    this_r ++;
	    next_r ++;
	    this_fr ++;
	    next_fr ++;
	    this_tag ++;
	    next_tag ++;
	    if (goon == 1){
		goon = 0;
	    }
	}
	else {
	    next_r = r.insert (next_r, mid_r);
	    next_fr = fr.insert (next_fr, mid_fr);
	    next_tag = tag.insert (next_tag, 1);
	}
    }
//    printf("tol=%e\tn_point=%d\n",tol,r.size());
    
    return integration;
}



static double lambda1;
static double lambda2;
static double indx;

inline double func10 (double r, double theta)
{
    double rr = r*r;
    double rcos = r*cos(theta);
    double rsin = r*sin(theta);
    double rrcos = rcos * rcos;
    double rrsin = rsin * rsin;

    return r * rrcos * pow(1 - rr, indx) * exp (lambda1 * rrcos + lambda2 * rrsin);
}

int main ()
{

  lambda1 = 90;
  lambda2 = 90;
  indx = 50;

  std::cout.precision (16);
  
  double m0;
  for (int time = 0; time < 1; time ++){
    m0 = int_1 (func10, 0, 1, 0, 2*M_PI, 1e-10, 20);
  }
  
  std::cout << "result is " << m0 << std::endl;
  
  return 0;
}


