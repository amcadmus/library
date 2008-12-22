#include "Stopwatch.h"
#include <iostream>

inline double Stopwatch::user () const
{
  return (double)(toc.tms_utime - tic.tms_utime) * HZi;
}

inline double Stopwatch::system () const
{
  return (double)(toc.tms_stime - tic.tms_stime) * HZi;
}

inline double Stopwatch::real () const
{
  return (double)(r1 - r0) * HZi;
}

inline void Stopwatch::stop ()
{
  r1 = times (&toc);
}

inline void Stopwatch::start() 
{
  r0 = times (&tic);
}

Stopwatch::Stopwatch ()
    : HZi (1./HZ)
{
}

int main(int argc, char * argv[])
{
  Stopwatch sw;
  double a = 1000;
  double b = 19;
  double c;
  sw.start();
  for (double i = 0; i < 1e8; i = i + 1.) c = a * b;
  sw.stop();  
  std::cout << "user time " << sw.user() << std::endl;
  std::cout << "system time " << sw.system() << std::endl;
  std::cout << "real time " << sw.real() << std::endl;

//   sw.start();
  for (double i = 0; i < 1e8; i = i + 1.) c = a * b;
  sw.stop();  
  std::cout << "user time " << sw.user() << std::endl;
  std::cout << "system time " << sw.system() << std::endl;
  std::cout << "real time " << sw.real() << std::endl;
  
}
