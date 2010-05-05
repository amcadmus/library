#include"Poly.h"
#include"Interpolation.h"
// #include"../../Software/include/ToolBox.h"
// #include"../../Software/src/ToolBox.h
#include <iostream>
#include <iterator>
#include <vector> 
#include <iostream>
#include <stdio.h>
#include <cmath>

int main(int argc, char * argv[])
{
//   std::vector<double > l, mu;
//   l.push_back(0.1);
//   l.push_back(0.2);
//   l.push_back(0.3);
//   l.push_back(0.4);
//   mu.push_back(1.2);
//   mu.push_back(3.4);
//   mu.push_back(5.6);
//   mu.push_back(7.8);
  

//   Interpolation::solverForSplinePeriodic (l.begin(), l.end(), mu.begin(), mu.end());
//   std::copy (mu.begin(), mu.end(), std::ostream_iterator<double > (std::cout, "\n"));

  // // FILE * fp = fopen ("averageHead.out", "r");
  // // double tmp1, tmp2;
  // // std::vector<double > x, y;
  // // while (fscanf(fp, "%lf%lf", &tmp1, &tmp2) != EOF){
  // //   x.push_back (tmp1);
  // //   y.push_back (tmp2);
  // // }
  // // x.push_back (x[2] - x[1] + x.back());
  // // y.push_back (y.front());
  // // PiecewisePoly ps;
  // // Interpolation::splinePeriodic (x, y, ps);
  
  // // for (unsigned i = 0; i < x.size()-1; ++i){
  // //   {
  // //     Poly tmp = ps.get_p()[i];
  // //     double v = tmp.value (x[i]);
  // //     tmp.derivative();
  // //     double dv = tmp.value (x[i]);
  // //     tmp.derivative();
  // //     double ddv = tmp.value (x[i]);
  // //     std::cout << v << '\t'
  // // 		<< dv << '\t'
  // // 		<< ddv << '\t' << std::endl;
  // //   }{
  // //     Poly tmp = ps.get_p()[i];
  // //     double v = tmp.value (x[i+1]);
  // //     tmp.derivative();
  // //     double dv = tmp.value (x[i+1]);
  // //     tmp.derivative();
  // //     double ddv = tmp.value (x[i+1]);
  // //     std::cout << v << '\t'
  // // 		<< dv << '\t'
  // // 		<< ddv << '\t' << std::endl;
  // //   }
  // // }
  

  std::vector<double > r;
  std::vector<double > ref;
  double rlo(0), rup(M_PI);
  unsigned Nr(100);
  double hr = (rup - rlo) / double(Nr);
  
  for (unsigned i = 0; i < Nr; ++i){
    r.push_back (rlo + i * hr);
    ref.push_back (sin(r.back()));
  }

  PiecewisePoly ps;
  Interpolation::piecewiseLinear (r, ref, ps);

  std::vector<double > x, y;
  double xlo(-4), xup(10);
  unsigned Nx (1000);
  double hx = (xup - xlo) / double(Nx);
  for (unsigned i = 0; i < Nx; ++i){
    x.push_back (xlo + i * hx);
  }
  ps.value_periodic (x, y);
  
  for (unsigned i = 0; i < Nx; ++i){
    std::cout << x[i] << " " << y[i] << std::endl;
  }
  for (unsigned i = 0; i < Nr; ++i){
    std::cerr << r[i] << " " << ref[i] << std::endl;
  }
  
  return 0;
}



// int main ()
// {
//   double low = 0;
//   double up = M_PI ;
//   int NN = 10;
//   double h = (up - low) / NN;
  
//   std::vector<double > x;
//   std::vector<double > y;
//   for (int i = 0; i < NN+1; ++i) {
//     double tmp = low;
//     if (i ==0 || i == NN) 
//       tmp += i * h;
//     else 
//       tmp += (i + (genrand_real1()-0.5)*0.9) * h;
//     x.push_back(tmp);
//     y.push_back(2*sin(tmp)+1);
//   }
  
//   PiecewisePoly p;
// //   Interpolation::spline (x, y, p);
//   Interpolation::spline (x.begin(), x.end(), y.begin(), p);
//   PiecewisePoly dp (p);
//   for (std::vector<Poly >::iterator i = dp.get_p().begin(); i != dp.get_p().end(); ++i){
//     i->derivative();
//   }
//   PiecewisePoly ddp (dp);
//   for (std::vector<Poly >::iterator i = ddp.get_p().begin(); i != ddp.get_p().end(); ++i){
//     i->derivative();
//   }

//   int Nr = 10 * N;
//   double hr = (up - low) / Nr;
//   std::vector<double > r;
//   std::vector<double > vr;
//   std::vector<double > dr;
//   std::vector<double > ddr;
//   for (int i = 0; i < Nr+1; ++i){
//     r.push_back( low + hr * i);
//   }
//   p.value (r, vr);
//   dp.value (r, dr);
//   ddp.value (r, ddr);
  
//   for (unsigned i = 0; i < r.size(); ++i){
//     std::cout << r[i] << '\t' 
// 	      << vr[i]  << '\t'
// 	      << dr[i]  << '\t'
// 	      << ddr[i]  << '\t'
// 	      << 1+2 * sin(r[i]) << '\t'
// 	      << 2 * cos(r[i]) << '\t'
// 	      << -2 * sin(r[i]) 
// 	      << std::endl;
//   }
  


//   return 0;
// }




// int main ()
// {
//   PiecewisePoly ps;
//   ps.clear();
//   ps.get_x().push_back(0);
//   ps.get_x().push_back(1);
//   ps.get_x().push_back(2);
//   ps.get_x().push_back(3);
//   Poly p;
//   p.one();
//   ps.get_p().push_back(p);
//   p *= 2;
//   ps.get_p().push_back(p);
//   p *= 3./2.;
//   ps.get_p().push_back(p);
  
// //   if (ps.valid())
// //     std::cout << "valid" << std::endl;
// //   else 
// //     std::cout << "invalid" << std::endl;
//   double r = 5.5;
//   std::vector<double > rs;
//   for (unsigned i = 0; i < 500; ++i){
//     rs.push_back(-1 + 5./500 * i);
//   }
//   std::vector<double > y;
//   ps.value (rs, y);

// //   std::cout << "single test: " << ps.value (r) << std::endl;
// //   std::cout << "multiple test: " << std::endl;
//   for (unsigned i = 0; i < y.size(); ++i){
//     std::cout << rs[i] << '\t' << y[i] << std::endl;
//   }
//   return 0;
// }

  
