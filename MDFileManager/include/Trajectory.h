#ifndef __MDFileManager_Trajectory_h_wanghan__
#define __MDFileManager_Trajectory_h_wanghan__

// #include "Defines.h"
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_xtc.h"
#include <vector> 

using namespace std;

class TrajLoader 
{
  XDRFILE *xd;
  int natoms;
  int step;
  float time;
  vector<double > box;
  rvec * xx;
  float prec;
  bool inited;
  void clear ();
public:
  TrajLoader ();
  ~TrajLoader ();
  TrajLoader (const char * filename);
  bool reinit (const char * filename);
  bool load ();
public:
  const vector<double > & getBox () const {return box;}
  float getTime () const {return time;}
public:
  void getFrame (vector<vector<double > > & frame);
}
    ;

#endif
