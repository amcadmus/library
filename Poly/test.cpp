#include"Poly.h"
#include"Interpolation.h"
#include"../../Software/include/ToolBox.h"
#include"../../Software/src/ToolBox.cpp"

int main ()
{
  double low = 0;
  double up = M_PI ;
  int NN = 10;
  double h = (up - low) / NN;
  
  std::vector<double > x;
  std::vector<double > y;
  for (int i = 0; i < NN+1; ++i) {
    double tmp = low;
    if (i ==0 || i == NN) 
      tmp += i * h;
    else 
      tmp += (i + (genrand_real1()-0.5)*0.9) * h;
    x.push_back(tmp);
    y.push_back(2*sin(tmp)+1);
  }
  
  PiecewisePoly p;
//   Interpolation::spline (x, y, p);
  Interpolation::spline (x.begin(), x.end(), y.begin(), p);
  PiecewisePoly dp (p);
  for (std::vector<Poly >::iterator i = dp.get_p().begin(); i != dp.get_p().end(); ++i){
    i->derivative();
  }
  PiecewisePoly ddp (dp);
  for (std::vector<Poly >::iterator i = ddp.get_p().begin(); i != ddp.get_p().end(); ++i){
    i->derivative();
  }

  int Nr = 10 * N;
  double hr = (up - low) / Nr;
  std::vector<double > r;
  std::vector<double > vr;
  std::vector<double > dr;
  std::vector<double > ddr;
  for (int i = 0; i < Nr+1; ++i){
    r.push_back( low + hr * i);
  }
  p.value (r, vr);
  dp.value (r, dr);
  ddp.value (r, ddr);
  
  for (unsigned i = 0; i < r.size(); ++i){
    std::cout << r[i] << '\t' 
	      << vr[i]  << '\t'
	      << dr[i]  << '\t'
	      << ddr[i]  << '\t'
	      << 1+2 * sin(r[i]) << '\t'
	      << 2 * cos(r[i]) << '\t'
	      << -2 * sin(r[i]) 
	      << std::endl;
  }
  


  return 0;
}




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

  
