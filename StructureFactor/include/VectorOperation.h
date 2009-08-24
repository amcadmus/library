#ifndef __VectorOperation_h_wanghan__
#define __VectorOperation_h_wanghan__

#include <vector>
#include <string>
#include <cmath>

namespace VectorOperation{
    
    inline double dot3 (const std::vector<double > & u, const std::vector<double > & v)
    {
      return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    }
    inline double dot2 (const std::vector<double > & u, const std::vector<double > & v)
    {
      return u[0] * v[0] + u[1] * v[1];
    }

}
;

#endif
