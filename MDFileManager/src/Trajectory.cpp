#include "Trajectory.h"
#include <stdlib.h>
#include <string.h>
#include <iostream>

void TrajLoader::
clear ()
{
  if (inited){
    free (xx);
    xdrfile_close (xd);
    time = 0.;
    inited = false;
  }
}

TrajLoader::
TrajLoader ()
    : time (0.),
      box (vector<double > (3, 0.)),
      inited (false)
{
}

TrajLoader::
~TrajLoader ()
{
  clear ();
}

TrajLoader::
TrajLoader (const char * filename)
    : time (0.),
      box (vector<double > (3, 0.)),
      inited (false)
{
  reinit (filename);
}  

bool TrajLoader::
reinit (const char * filename)
{
  char tmpname[2048];
  strncpy (tmpname, filename, 2047);
  
  xd = xdrfile_open (filename, "r");
  if (xd == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    return false;
  }
  read_xtc_natoms (tmpname, &natoms);
  step = 0;
  time = 0.;
  box = vector<double > (3, 0.);
  xx = (rvec *) malloc (sizeof(rvec) * natoms);
  prec = 1000.;

  inited = true;

  load ();

  xdrfile_close (xd);
  xd = xdrfile_open (filename, "r");  
  if (xd == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    clear ();
    return false;
  }
  
  return true;
}

bool TrajLoader::
load ()
{
  if (inited){
    matrix tmpBox;
    int st = read_xtc (xd, natoms, &step, &time, tmpBox, xx, &prec);
    box[0] = tmpBox[0][0];
    box[1] = tmpBox[1][1];
    box[2] = tmpBox[2][2];
    if (st == exdrOK) return true;
    else return false;
  }
  else {
    std::cerr << "not initiated, do nothing." << std::endl;
    return false;
  }
}

void TrajLoader::
getFrame (vector<vector<double > > & frame)
{
  frame.resize (natoms);
  for (int ii = 0; ii < natoms; ++ii){
    frame[ii].resize(3);
    frame[ii][0] = xx[ii][0];
    frame[ii][1] = xx[ii][1];
    frame[ii][2] = xx[ii][2];
  }
}


