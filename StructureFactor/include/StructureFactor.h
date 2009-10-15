#ifndef __StructureFacture_h_wanghan__
#define __StructureFacture_h_wanghan__

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <queue>
#include <vector>
#include <cmath>
#include "InterpolationBases.h"
#include <fftw3.h>
#include <complex>

typedef double value_type;

template <int DomainDim>
class StructureFactor 
{
  std::vector<unsigned > K;
  std::vector<std::vector<value_type > > vecA;
  std::vector<std::vector<value_type > > vecAStar;
  smoothInterpolationBase * Mn;
  value_type V;
  std::vector<std::vector<std::complex<double > > > valueb;
  void build ();
private:
  fftw_complex * Q;	// do not forget to set the image part to 0
  fftw_complex * QF;
  fftw_plan forwardQ;
private:
  bool check_parameter_consistency ();
  void calV ();
  void calAStar ();
  void calB ();
  void calQ (const std::vector<std::vector<double > > & coord,
	     const std::vector<double > & value);
public:
  ~StructureFactor();
  void init (const std::vector<std::vector<double > > &vecA_,  
// vecA is the box base vector
	     const std::vector<unsigned > K_,
	     const InterpolationBaseInfo::Order & interpOrder_);

  void calStructureFactor(
      const std::vector<std::vector<double > > & coord,
      const std::vector<double > & value,
      std::vector<std::complex<double > > & sf);
  void test(
      const std::vector<std::vector<double > > & coord,
      const std::vector<double > & value,
      std::vector<std::complex<double > > & sf);
} 
    ;

namespace StructureFactorUtils {
inline double dot3 (const std::vector<double > & u, const std::vector<double > & v)
{
  return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}
inline double dot2 (const std::vector<double > & u, const std::vector<double > & v)
{
  return u[0] * v[0] + u[1] * v[1];
}
};


template <int DomainDim>
void StructureFactor<DomainDim >::calQ (
    const std::vector<std::vector<double > > & coord,
    const std::vector<double > & value)
{
  //   std::cerr << "K0  " << K[0] << std::endl;
  unsigned n = Mn->getN();

  if (DomainDim == 3){
    for (unsigned i = 0; i < K[0] * K[1] * K[2]; i ++){
      Q[i][1] = 0;
      Q[i][0] = 0;
    }
    double ii0 = 1./ value_type(K[0]);
    double ii1 = 1./ value_type(K[1]);
    double ii2 = 1./ value_type(K[2]);
    for (unsigned i = 0; i < coord.size(); ++i){
      std::vector<value_type > u(3);
      u[0] = K[0] * StructureFactorUtils::dot3 (vecAStar[0], coord[i]);
      u[1] = K[1] * StructureFactorUtils::dot3 (vecAStar[1], coord[i]);
      u[2] = K[2] * StructureFactorUtils::dot3 (vecAStar[2], coord[i]);
      int A0 = -int(floor ((u[0]) * ii0)) ;
      int A1 = -int(floor ((u[1]) * ii1)) ;
      int A2 = -int(floor ((u[2]) * ii2)) ;
      value_type posi0 = u[0] + A0 * K[0];
      value_type posi1 = u[1] + A1 * K[1];
      value_type posi2 = u[2] + A2 * K[2];
      value_type tmp0 = 0;
      value_type tmp1 = 0;
      value_type tmp2 = 0;

      unsigned index0, index1;
      for (int k0 = int(ceil(posi0-n)); k0 < int(ceil(posi0)); ++ k0){
	Mn->value (posi0 - k0, tmp0);
	index0 = K[1] * (k0<0 ? k0+K[0] : k0);
	for (int k1 = int(ceil(posi1-n)); k1 < int(ceil(posi1)); ++ k1){
	  Mn->value (posi1 - k1, tmp1);
	  index1 = K[2] * ((k1<0 ? k1+K[1] : k1) + index0);
	  for (int k2 = int(ceil(posi2-n)); k2 < int(ceil(posi2)); ++ k2){
	    Mn->value (posi2 - k2, tmp2);
	    Q[(k2<0 ? k2+K[2] : k2) + index1][0] += value[i] * tmp0 * tmp1 * tmp2;
	  }
	}
      }
    }
  }
  else if (DomainDim == 2) {
    for (unsigned i = 0; i < K[0] * K[1]; i ++){
      Q[i][1] = 0;
      Q[i][0] = 0;
    }
    double ii0 = 1./ value_type(K[0]);
    double ii1 = 1./ value_type(K[1]);
    for (unsigned i = 0; i < coord.size(); ++i){
      std::vector<value_type > u(2);
      u[0] = K[0] * StructureFactorUtils::dot2 (vecAStar[0], coord[i]);
      u[1] = K[1] * StructureFactorUtils::dot2 (vecAStar[1], coord[i]);
      int A0 = -int(floor ((u[0]) * ii0)) ;
      int A1 = -int(floor ((u[1]) * ii1)) ;
      value_type posi0 = u[0] + A0 * K[0];
      value_type posi1 = u[1] + A1 * K[1];
      value_type tmp0 = 0;
      value_type tmp1 = 0;

      unsigned index0;
      for (int k0 = int(ceil(posi0-n)); k0 < int(ceil(posi0)); ++ k0){
	Mn->value (posi0 - k0, tmp0);
	index0 = K[1] * (k0<0 ? k0+K[0] : k0);
	for (int k1 = int(ceil(posi1-n)); k1 < int(ceil(posi1)); ++ k1){
	  Mn->value (posi1 - k1, tmp1);
	  Q[((k1<0 ? k1+K[1] : k1) + index0)][0] += value[i] * tmp0 * tmp1;
	}
      }
    }
  }
}

template <int DomainDim>
void StructureFactor<DomainDim >::calStructureFactor(
    const std::vector<std::vector<double > > & coord,
    const std::vector<double > & value,
    std::vector<std::complex<double > > & sf)
{
  calQ (coord, value);
  fftw_execute (forwardQ);
  sf.clear();
  
  if (DomainDim == 3){
    for (unsigned m0 = 0; m0 < K[0]; ++m0){
      for (unsigned m1 = 0; m1 < K[1]; ++m1){
	for (unsigned m2 = 0; m2 < K[2]; ++m2){
	  std::complex<double > tmp (QF[m2 + K[2] * (m1 + K[1] * m0)][0],
				     QF[m2 + K[2] * (m1 + K[1] * m0)][1]);
	  sf.push_back (valueb[0][m0] * valueb[1][m1] * valueb[2][m2] * tmp);
	}
      }
    }
  }
  else if (DomainDim == 2){
    for (unsigned m0 = 0; m0 < K[0]; ++m0){
      for (unsigned m1 = 0; m1 < K[1]; ++m1){
	std::complex<double > tmp (QF[ (m1 + K[1] * m0)][0],
				   QF[ (m1 + K[1] * m0)][1]);
	sf.push_back (valueb[0][m0] * valueb[1][m1]  * tmp);
      }
    }
  }
}

template <int DomainDim>
void StructureFactor<DomainDim>::init (
    const std::vector<std::vector<value_type > > &vecA_,
    const std::vector<unsigned > K_,
    const InterpolationBaseInfo::Order & interpOrder_)
{
  switch (interpOrder_){
  case 2 :
      Mn = new BSpline2();
      break;
  case 4 :
      Mn = new BSpline4();
      break;
  case 6 :
      Mn = new BSpline6();
      break;
  case 8 :
      Mn = new BSpline8();
      break;
  case 10 :
      Mn = new BSpline10();
      break;
  case 12 :
      Mn = new BSpline12();
      break;
  case 14 :
      Mn = new BSpline14();
      break;
  default:
      std::cerr << "no such order implemented, use order 8" << std::endl;
      Mn = new BSpline8();
  }
  K = K_;
  vecA = vecA_;
  bool status_pc = check_parameter_consistency ();
  if (!status_pc ){
    std::cerr << "StructureFactor: one of the initalize parameter is invalid!\n" ;
    exit (1);
  }
  calV();
  calAStar();
  calB ();

  build ();  
}

template <int DomainDim >
void StructureFactor<DomainDim >::build()
{

  if (DomainDim == 3){
    int size = K[0] * K[1] * K[2];
    Q	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);
    QF	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);

    forwardQ	= fftw_plan_dft_3d (K[0], K[1], K[2], Q  , QF, 1  , FFTW_MEASURE);
  }
  else if (DomainDim == 2){
    int size = K[0] * K[1];
    Q	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);
    QF	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);

    forwardQ	= fftw_plan_dft_2d (K[0], K[1], Q  , QF, 1  , FFTW_MEASURE);
  }  
}


template <int DomainDim>
void StructureFactor<DomainDim >::calB ()
{
  unsigned n = Mn->getN();

  valueb.clear();
  valueb.resize (DomainDim);
  for (unsigned i = 0; i < DomainDim; ++i){
    valueb[i].resize (K[i]);
  }
  for (unsigned i = 0; i < DomainDim; ++i){
    for (unsigned m = 0; m < K[i]; ++m){
      std::vector<value_type > fenzi (2);
      value_type tmp = 2 * M_PI * (n-1) * m / value_type (K[i]);
      fenzi[0] = cos(tmp);
      fenzi[1] = sin(tmp);
      std::vector<value_type > fenmu (2, 0);
      for (unsigned k = 0; k < n-1; k ++){
	value_type scale ;
	Mn->value (k+1, scale);
	tmp = 2 * M_PI * m * k / value_type (K[i]);
	fenmu[0] += scale * cos(tmp);
	fenmu[1] += scale * sin(tmp);
      }
      std::vector<value_type > btmp (2);
      value_type scale = 1./ (fenmu[0]*fenmu[0] + fenmu[1]*fenmu[1]);
      btmp[0] = scale * (fenzi[0] * fenmu[0] + fenzi[1] * fenmu[1]);
      btmp[1] = scale * (fenzi[1] * fenmu[0] - fenzi[0] * fenmu[1]);
      valueb[i][m] = std::complex<double >(btmp[0], btmp[1]);
    }
  }  
}



template <int DomainDim>
bool StructureFactor<DomainDim>::check_parameter_consistency ()
{
  if (DomainDim == 2){
    if (K.size() != 2) return false;
    if (vecA.size() != 2) return false;
    if (vecA[0].size() != 2) return false;
    if (vecA[1].size() != 2) return false;
  }
  else if (DomainDim == 3){
    if (K.size() != 3) return false;
    if (vecA.size() != 3) return false;
    if (vecA[0].size() != 3) return false;
    if (vecA[1].size() != 3) return false;
    if (vecA[2].size() != 3) return false;
  }
  else {
    return false;
  }
  int n = Mn->getN();
  for (int i = 0; i < DomainDim; ++i){
    if (unsigned (n) >= K[i]){
      std::cerr << "StructureFactor: interpolation order should be biger than number of grid points" << std::endl;
      return false;
    }
  }

  return true;
}

template <int DomainDim>
void StructureFactor<DomainDim>::calV()
{
  if (DomainDim == 3){
    V = vecA[0][0] * (vecA[1][1]*vecA[2][2] - vecA[2][1]*vecA[1][2]) - 
	vecA[0][1] * (vecA[1][0]*vecA[2][2] - vecA[2][0]*vecA[1][2]) +
	vecA[0][2] * (vecA[1][0]*vecA[2][1] - vecA[2][0]*vecA[1][1]);
  }
}

template <int DomainDim>
void StructureFactor<DomainDim>::calAStar ()
{
  if (DomainDim == 3){
    vecAStar.resize (3);
    vecAStar[0].resize (3);
    vecAStar[1].resize (3);
    vecAStar[2].resize (3);
    vecAStar[0][0] =( vecA[1][1]*vecA[2][2] - vecA[2][1]*vecA[1][2]) / V;
    vecAStar[1][1] =( vecA[0][0]*vecA[2][2] - vecA[2][0]*vecA[0][2]) / V;
    vecAStar[2][2] =( vecA[0][0]*vecA[1][1] - vecA[1][0]*vecA[0][1]) / V;
    vecAStar[1][0] =(-vecA[1][0]*vecA[2][2] + vecA[2][0]*vecA[1][2]) / V;
    vecAStar[2][0] =( vecA[1][0]*vecA[2][1] - vecA[2][0]*vecA[1][1]) / V;
    vecAStar[0][1] =(-vecA[0][1]*vecA[2][2] + vecA[2][1]*vecA[0][2]) / V;
    vecAStar[2][1] =(-vecA[0][0]*vecA[2][1] + vecA[2][0]*vecA[0][1]) / V;
    vecAStar[0][2] =( vecA[0][1]*vecA[1][2] - vecA[1][1]*vecA[0][2]) / V;
    vecAStar[1][2] =(-vecA[0][0]*vecA[1][2] + vecA[1][0]*vecA[0][2]) / V;
  }
  else if (DomainDim == 2){
    vecAStar.resize (2);
    vecAStar[0].resize (2);
    vecAStar[1].resize (2);
    std::vector<double > a0p (2), a1p(2);
    a0p[0] =  vecA[0][1];
    a0p[1] = -vecA[0][0];
    a1p[0] =  vecA[1][1];
    a1p[1] = -vecA[1][0];
    double scale1 = StructureFactorUtils::dot2 (a0p, vecA[1]);
    double scale0 = StructureFactorUtils::dot2 (a1p, vecA[0]);
    vecAStar[0][0] = a1p[0] / scale0;
    vecAStar[0][1] = a1p[1] / scale0;
    vecAStar[1][0] = a0p[0] / scale1;
    vecAStar[1][1] = a0p[1] / scale1;
    std::cout << StructureFactorUtils::dot2 (vecAStar[0], vecA[0]) << std::endl;
    std::cout << StructureFactorUtils::dot2 (vecAStar[0], vecA[1]) << std::endl;
    std::cout << StructureFactorUtils::dot2 (vecAStar[1], vecA[0]) << std::endl;
    std::cout << StructureFactorUtils::dot2 (vecAStar[1], vecA[1]) << std::endl;
    std::cout << vecAStar[0][0] << " " << vecAStar[0][1] << std::endl;
    std::cout << vecAStar[1][0] << " " << vecAStar[1][1] << std::endl;
    
  }
}

template <int DomainDim > 
StructureFactor<DomainDim >::~StructureFactor()
{
  delete Mn;
  fftw_free (Q);
  fftw_free (QF);
  fftw_destroy_plan (forwardQ);
}

template <int DomainDim>
void StructureFactor<DomainDim>::test(
    const std::vector<std::vector<double > > & coord,
    const std::vector<double > & value,
    std::vector<std::complex<double > > & sf)
{
  fftw_complex * in, *ot;
  fftw_plan forward;
  
  int size = K[0] * K[1] * K[2];
  in	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);
  ot	= (fftw_complex *) fftw_malloc (sizeof(fftw_complex) * size);
  
  forward = fftw_plan_dft_3d (K[0], K[1], K[2], in  , ot, 1  , FFTW_MEASURE);

  for (unsigned i = 0; i < value.size(); ++i){
    in[i][0] = value[i];
    in[i][1] = 0;
  }
  
  fftw_execute (forward);

  sf.clear();
  for (unsigned i = 0; i < value.size(); ++i){
    std::complex<double > tmp (ot[i][0], ot[i][1]);
    sf.push_back(tmp);
  }

  fftw_free (in);
  fftw_free (ot);
  fftw_destroy_plan (forward);
  
}




#endif
