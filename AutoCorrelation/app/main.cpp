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

#include "PieceTraj.h"
#include "AutoCorrelCalculator.h"
#include "Trajectory.h"

int main(int argc, char * argv[])
{
  // PieceTraj<double > traj;
  // traj.push_back(1.);

  TrrLoader trr ("traj.trr");

  int natoms = trr.getNatoms();
  int nFrame = 100;

  vector<AutoCorrelCalculator> accx (natoms);
  vector<AutoCorrelCalculator> accy (natoms);
  vector<AutoCorrelCalculator> accz (natoms);
  vector<double > times;
  BlockAverage_acc avgVar;
  double sum = 0.;
  int countFrame = 0;
  for (int ii = 0; ii < natoms; ++ii){
    accx[ii].reinit (nFrame, 100);
    accy[ii].reinit (nFrame, 100);
    accz[ii].reinit (nFrame, 100);
  }

  vector<vector<double > > xx, vv, ff;  
  while (trr.load ()){    
    if (trr.getTime() < 10) continue;
    times.push_back (trr.getTime());
    // printf ("# read at time %f\n", trr.getTime());
    // if (trr.getTime() > 20) break;
    trr.getFrame (xx, vv, ff);
    for (unsigned ii = 0; ii < vv.size(); ++ii){
      // water, skip atom hydrogen
      if (ii % 3 != 0) continue;
      accx[ii].push_back (vv[ii][0]);
      accy[ii].push_back (vv[ii][1]);
      accz[ii].push_back (vv[ii][2]);
      avgVar.deposite (vv[ii][0] * vv[ii][0] + vv[ii][1] * vv[ii][1] + vv[ii][2] * vv[ii][2]);
      sum += vv[ii][0] * vv[ii][0] + vv[ii][1] * vv[ii][1] + vv[ii][2] * vv[ii][2];
    }
    countFrame ++;
  }

  for (int ii = 0; ii < natoms; ++ii){
      // water, skip atom hydrogen
    if (ii % 3 != 0) continue;
    accx[ii].calculate ();
    accy[ii].calculate ();
    accz[ii].calculate ();
  }
  avgVar.calculate();
  printf ("# var vel: %f %f\n", avgVar.getAvg(), sum / natoms / countFrame);

  vector<BlockAverage_acc > correl (accx[0].nData());
  
  for (unsigned ii = 0; ii < correl.size(); ++ii){
    for (int jj = 0; jj < natoms; ++jj){
      // water, skip atom hydrogen
      if (jj % 3 != 0) continue;
      correl[ii].deposite (accx[jj].value(ii));
      correl[ii].deposite (accy[jj].value(ii));
      correl[ii].deposite (accz[jj].value(ii));
    }
  }

  for (unsigned ii = 0; ii < correl.size(); ++ii){
    correl[ii].calculate();
    printf ("%f %e %e\n",
	    times[ii]-times[0],
	    correl[ii].getAvg() * 3,
	    correl[ii].getAvgError() * 3);
  }
}
