#ifndef __Target_h_wanghan__
#define __Target_h_wanghan__

typedef  double ScalorType;

class OPTIMIZATION_TARGET 
{
  ScalorType hg;
  ScalorType hh;
public:
  OPTIMIZATION_TARGET ()
      : hg (1e-6), hh(1e-4) {}
  void value (const std::vector<ScalorType > & x,
	      ScalorType & value) const;
  void grad  (const std::vector<ScalorType > & x,
	      std::vector<ScalorType > & g) const;
  void hessian (const std::vector<ScalorType > & x,
		std::vector<ScalorType > & hessian) const;
};

#endif
