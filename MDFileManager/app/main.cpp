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

#include "Trajectory.h"

using namespace std;

int main(int argc, char * argv[])
{
  TrajLoader traj;

  std::cout << traj.getBox()[0] << std::endl;
  std::cout << traj.getBox()[1] << std::endl;
  std::cout << traj.getBox()[2] << std::endl;
  std::cout << traj.getTime() << std::endl;

  traj.reinit ("traj.xtc");

  std::cout << endl;
  std::cout << traj.getBox()[0] << std::endl;
  std::cout << traj.getBox()[1] << std::endl;
  std::cout << traj.getBox()[2] << std::endl;
  std::cout << traj.getTime() << std::endl;

  traj.load ();

  std::cout << endl;
  std::cout << traj.getBox()[0] << std::endl;
  std::cout << traj.getBox()[1] << std::endl;
  std::cout << traj.getBox()[2] << std::endl;
  std::cout << traj.getTime() << std::endl;

  // traj.load ();
  // traj.load ();
  
  vector<vector<double > > frame;
  traj.getFrame (frame);
  for (unsigned ii = 0; ii < frame.size(); ++ii){
    printf ("%d: %f   %f %f %f\n", ii, traj.getTime(), frame[ii][0], frame[ii][1], frame[ii][2]);
  }
  
  return 0;
}
