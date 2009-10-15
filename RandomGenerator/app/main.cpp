#include "RandomGenerator.h"
#include<iostream>

int main(int argc, char * argv[])
{
  double a = RandomGenerator_MT19937::genrand_real1();
  std::cout << a << std::endl;
  return 0;
}
