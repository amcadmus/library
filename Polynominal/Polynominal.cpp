#include "Polynominal.h"

Term::Term (const std::string & name_, 
	    const double & scale_, 
	    const std::vector<int > & index_)
    : name (name_), scale (scale_), index (index_) 
{
}

Term::Term (const std::string & name_, const unsigned & dim, const double & scale_)
    : name (name_), scale (scale_), index (dim, 0)
{
}

Term::Term (const std::string & name_, const unsigned & dim)
    : name (name_), scale (0.0), index (dim, 0)
{
}

bool Term::isZero ()
{
  return (scale == 0.);
}

void Term::regZero ()
{
  if (isZero()) 
    std::fill (index.begin(), index.end(), 0);
}

void Term::reinit (const std::string & name_, 
		   const double & scale_, 
		   const std::vector<int > & index_)
{
  name = name_;
  scale = scale_;
  index = index_;
}
  
Term & Term::operator += (const Term & term)
{
  if (isZero()){
    *this = term;
  }
  else if (term.name == name){
    if (term.index == index){
      scale += term.scale;
    }
  }
  return *this;
}

Term & Term::operator *= (const Term & term)
{
  if (term.name == name){
    scale *= term.scale;
    std::vector<int >::iterator i = index.begin();
    std::vector<int >::const_iterator j = term.index.begin();
    for (; i != index.end(); ++i, ++j){
      *i += *j;
    }
  }
  return *this;
}

Term & Term::operator *= (const double & scale_)
{
  scale *= scale;
  return *this;
}

void Term::derivative (const unsigned & i)
{
  scale *= (index[i] --);
  if (index[1] == -1){
    index[1] = 0;
  }
}


void Term::print ()
{
  std::cout << scale << " * " ;
  for (unsigned i = 0; i < index.size(); i ++){
    std::cout << name << "_" << i << "^" << index[i];
  }
}

Polynominal::Polynominal (const std::string & name_)
    : name (name_)
{
}

Polynominal & Polynominal::operator += (const Term & term)
{
  if (term.getName () != name){
    std::cout << 0. << std::endl;
    return *this;
  }
  if (this->empty()){
    this->push_back (term);
    return * this;
  }
  for (std::list<Term >::iterator i = this->begin(); i != this->end(); ++ i){
    if ( (i->getIndex() < term.getIndex())){
      if (i->getIndex() == term.getIndex()){
	i->operator += (term);
	return * this;
      }
      else {
	this->insert (i, term);
	return *this;
      }
    }
  }
  this->push_back (term);
  return *this;
}

Polynominal & Polynominal::operator *= (const Term & term)
{
  for (std::list<Term >::iterator i = this->begin(); i != this->end(); ++ i){
    (*i) *= term;
  }
  return *this;
}
 
void Polynominal::print ()  
{
  if (! this->std::list<Term >::empty()){
    std::list<Term >::iterator i = this->begin();
    i->print();
    ++ i;
    for (; i != this->end(); ++ i){
      std::cout << " + ";
      i->print();
    }
    std::cout << std::endl;
  }
}

void Polynominal::derivative (const unsigned & i)
{
  for (std::list<Term >::iterator it = this->begin(); it != this->end(); ++ it){
    it->getScale() *= (it->getIndex(i) --);
    if (it->getIndex(i) == -1){
      it->getIndex(i) = 0;
    }
  }
  delZero();
}

void Polynominal::delZero ()
{
  for (std::list<Term >::iterator i = this->begin(); i != this->end();){
    if (i->isZero()){
      i = this->erase (i);
      if (i == this->end()){
	return;
      }
    }
    else{
      ++ i;
    }
  }
}

// int main(int argc, char * argv[])
// {
//   std::vector<int > index (3,0);
//   Polynominal poly ("Q");
//   Term tmp ("Q", 3);
//   Term tmp1 ("Q", 3);

//   std::cout << index.size() << std::endl;

//   index[0] = 3;
//   tmp.reinit ("Q", 05, index);
//   tmp.print ();
//   std::cout << std::endl;
//   poly += (tmp);
//   poly.print ();

//   index[0] = 1;
//   tmp.reinit ("Q", 1.5, index);
//   poly += (tmp);
//   poly.print ();
  
//   index[0] = 2;
//   index[1] = 1;
//   tmp.reinit ("Q", 3.5, index);
//   poly += (tmp);
//   poly.print ();
  
//   poly.delZero();
//   poly.print();

//   poly.derivative (1);
//   poly.print();

//   tmp.reinit ("Q", 3.5, index);
//   tmp.print();
//   std::cout << std::endl;
// //   index[2] = 2;
// //   tmp1.reinit ("Q", 2.5, index);
// //   tmp1.print();
// //   std::cout << std::endl;
// //   tmp1 *= tmp;
// //   tmp1.print();
// //   std::cout << std::endl;

//   poly *= tmp;
//   poly.print();

//   return 0;
// }







