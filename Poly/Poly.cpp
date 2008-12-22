#include "Poly.h"

bool PiecewisePoly::valid () const
{
  if (x.size() != p.size()+1) return false;
  std::vector<double >::const_iterator i = x.begin();
  std::vector<double >::const_iterator j = x.begin();
  for (++j ;j != x.end(); ++i, ++j){
    if (*i > *j) return false;
  }
  return true;
}

double PiecewisePoly::value (const double & xx) const
{
  unsigned begin = 0;
  unsigned end = x.size() - 1;
  unsigned mid = end/2;
  if (end == begin){
    return 0;
  }
  while (end - begin > 1){
    if (xx < x[mid]){
      end = mid;
      mid = (begin + end) / 2;
    }
    else{
      begin = mid;
      mid = (begin + end) / 2;
    }
  }
  return p[begin].value(xx);
}

double PiecewisePoly::value (const double & xx,
			     unsigned & begin,
			     unsigned & end) const
{
  if (end <= begin) return 0;
  if (end - begin == 1)  return p[begin].value(xx);
  unsigned mid = (begin + end)/2;
  while (end - begin > 1){
    if (xx < x[mid]){
      end = mid;
      mid = (begin + end) / 2;
    }
    else{
      begin = mid;
      mid = (begin + end) / 2;
    }
  }
  return p[begin].value(xx);
}

void PiecewisePoly::value (const unsigned & xbegin,
			   const unsigned & xend,
			   const std::vector<double > & r,
			   const unsigned & rbegin,
			   const unsigned & rend,
			   std::vector<double > & y) const
{
  unsigned xbegin1 = xbegin;
  unsigned xend1 = xend;
  if (rend - rbegin <= 1){
    y[rbegin] = value(r[rbegin], xbegin1, xend1);
    xbegin1 = xbegin;
    xend1 = xend;
    y[rend] = value(r[rend], xbegin1, xend1);
  }
  else {
    unsigned rmid = (rbegin + rend) / 2;
    y[rmid] = value (r[rmid], xbegin1, xend1);
    value (xbegin, xend1, r, rbegin, rmid-1, y);
    value (xbegin1, xend, r, rmid+1, rend, y);
  }
}

// suppose that
void PiecewisePoly::value (const std::vector<double > & r,
			   std::vector<double > & y) const
{
  y.resize(r.size());
  value (0, x.size()-1, r, 0, r.size()-1, y);
}

double Poly::value (const double & x) const
{
  double value = a[a.size()-1];
  for (int i = a.size() - 2; i >= 0; --i){
    value = value * x + a[i];
  }
  return value;
}

Poly::Poly ()
    : a (1, 0.)
{
  order = 0;
}

Poly::Poly (const std::vector<double > & out)
    : a(out) 
{
  order = out.size() - 1;
}

void Poly::reinit (const std::vector<double > & out)
{
  a = out;
  order = out.size() - 1;
} 

Poly & Poly::operator = (const Poly & p)
{
  a = p.a;
  order = p.order;
  return *this;
}

Poly & Poly::operator += (const Poly & p)
{
  if (p.a.size() > a.size()){
    a.resize(p.a.size(), 0);
    order = p.order;
    for (unsigned i = 0; i <= order; i ++){
      a[i] += p.a[i];
    }
  }
  else {
    for (unsigned i = 0; i <= p.order; i ++){
      a[i] += p.a[i];
    }
  }
  return * this;
}

Poly & Poly::derivative ()
{
  if (order == 0) {
    a[0] = 0;
    return *this;
  }
  for (unsigned i = 0; i < order; i ++){
    a[i] = a[i+1] * (i+1);
  }
  order --;
  a.pop_back();
  return * this;
}


Poly & Poly::operator *= (const double & scale)
{
  if (scale == 0){
    order = 0;
    a.resize (1);
    a[0] = 0;
  }
  else {
    for (std::vector<double >::iterator i = a.begin(); i != a.end(); i ++){
      *i *= scale;
    }
  }
  return * this;
}

  

Poly & Poly::operator *= (const Poly & p)
{
  std::vector<double > a1 (a);
  unsigned order1 (order);
  
  order += p.order;
  a.resize (order+1, 0);
  
  for (std::vector<double >::iterator i = a.begin(); i != a.end(); i ++){
    *i *= p.a[0];
  }
  if (p.order >= 1){
    for (unsigned i = 1; i <= p.order; i ++){
      for (unsigned j = 0; j <= order1; j ++){
	a[i+j] += a1[j] * p.a[i];
      }
    }
  }
  return *this;
}

void Poly::print ()
{
  for (unsigned i = 0; i <= order; i ++){
    std::cout << a[i] <<'\t' ;
  }
  std::cout << std::endl;
}

void Poly::print (const std::string & x)
{
  std::cout << a[0];
  for (unsigned i = 1; i <= order; i ++){
    std::cout << " + " <<  a[i] << x << "^" << i ;
  }
  std::cout << std::endl;
}

void Poly::printCode (const std::string & x)
{
  std::cout.precision (16);
  if (order == 0){
    std::cout << a[0] << std::endl;
    return;
  }
  
  for (unsigned i = 0; i < order-1; i ++){
    std::cout << "(" ;
  }
  std::vector<double >::reverse_iterator p = a.rbegin();
  std::cout << *(p++) << " * " << x << " + " ;
  std::cout << *(p++);
  for (; p != a.rend(); p ++){
    std::cout << ") * " << x << " + " << *p;
  }
  std::cout << std::endl;
}
