#ifndef __IntegerDecomposer_h_wanghan__
#define __IntegerDecomposer_h_wanghan__

#include<vector>

class IntegerDecomposer // find a integer that is very near to the
			// input and can be decomposed to 2,3,5 and
			// 7. This is very useful to choose a grid for
			// FFTW.
{
  std::vector<std::vector<int > > good; 
public:
  IntegerDecomposer ();
  int decompose (int input, std::vector<int > & index);
}
    ;

#endif
