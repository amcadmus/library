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
#include "StructureFactor.h"
#include "ToolBox.h"

void generateSystem (
    const std::vector<std::vector<double > > & vecA,
    const unsigned & N,
    std::vector<std::vector<double > > & coord,
    std::vector<double > & value)
{
  std::vector<double > box(3, 0);
  box[0] = vecA[0][0];
  box[1] = vecA[1][1];
  box[2] = vecA[2][2];
  coord.clear();
  value.clear();
  for (unsigned i = 0; i < N; ++i){
    std::vector<double > tmp (3, 0);
    tmp[0] = ToolBox::genrand_real2() * box[0];
    tmp[1] = ToolBox::genrand_real2() * box[1];
    tmp[2] = ToolBox::genrand_real2() * box[2];
    coord.push_back (tmp);
    // value.push_back (ToolBox::genrand_real2());
    // value.push_back (cos(tmp[0] / box[0] * 2 * M_PI));
    value.push_back(1);
  }
}

void generateSystem (
    const std::vector<std::vector<double > > & vecA,
    const std::vector<unsigned > K,
    std::vector<std::vector<double > > &coord,
    std::vector<double > & value )
{
  std::vector<double > box(3, 0);
  box[0] = vecA[0][0];
  box[1] = vecA[1][1];
  box[2] = vecA[2][2];
  coord.clear();
  value.clear();
  for (unsigned m0 = 0; m0 < K[0]; ++m0){
    for (unsigned m1 = 0; m1 < K[1]; ++m1){
      for (unsigned m2 = 0; m2 < K[2]; ++m2){
	std::vector<double > tmp(3, 0);
	tmp[0] = m0 / double(K[0]) * box[0];
	tmp[1] = m1 / double(K[1]) * box[1];
	tmp[2] = m2 / double(K[2]) * box[2];
	coord.push_back(tmp);
	// value.push_back(ToolBox::genrand_real2());
	value.push_back (cos(tmp[0] / box[0] * 2 * M_PI));
      }
    }
  }
}    

void test_sf (
    const std::vector<std::vector<double > > &vecA,
    const std::vector<int > m,
    const std::vector<std::vector<double > > & coord,
    const std::vector<double > & value,
    std::complex<double > & result)
{
  
  std::vector<double > box(3, 0);
  box[0] = vecA[0][0];
  box[1] = vecA[1][1];
  box[2] = vecA[2][2];

  result = std::complex<double > (0, 0);
  
  for (unsigned i = 0; i < coord.size(); ++i){
    double tmp = 2 *M_PI * (m[0]*coord[i][0] / box[0] +
			    m[1]*coord[i][1] / box[1] +
			    m[2]*coord[i][2] / box[2] );
			    
    result.real() += value[i] * cos(tmp);
    result.imag() += value[i] * sin(tmp);
  }
}


int main(int argc, char * argv[])
{
  StructureFactor <3> sf;
  std::vector<std::vector<double > > vecA (3);
  vecA[0].resize (3, 0);
  vecA[1].resize (3, 0);
  vecA[2].resize (3, 0);
  vecA[0][0] = 3;
  vecA[1][1] = 3;
  vecA[2][2] = 3;
  std::vector<unsigned > K(3, 32);
  // K[1] = 32;
  // K[2] = 32;

  
  std::vector<std::vector<double > > coord;
  std::vector<double > value;
  
  generateSystem (vecA, 100, coord, value);
  // std::cout << value.size() << std::endl;
  std::complex<double > result0 ;
  std::vector<int > testindx (3, 0);
  // testindx[0] = 10;
  testindx[1] = 7;
  testindx[2] = 7;
  testindx[0] = 7;
  
  test_sf (vecA, testindx, coord, value, result0);
  

  sf.init (vecA, K, 10);
  std::vector<std::complex<double > > result1;
  sf.calStructureFactor (coord, value, result1);
  
  std::cout << result0 << std::endl;
  std::cout << result1[testindx[2] + (K[2]) * (testindx[1] + K[1] * testindx[0])] << std::endl;


  // std::cout << value.size() << std::endl;
  
  // sf.test (coord, value, result1);
  // std::cout << result1[testindx[2] + (K[2]) * (testindx[1] + K[1] * testindx[0])] << std::endl;
  
  return 0;
  
}
