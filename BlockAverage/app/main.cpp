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

#include "BlockAverage.h"
using namespace std;

int main(int argc, char * argv[])
{
  FILE * fp;
  if ((fp = fopen("energy.xvg", "r")) == NULL) {
      std::cerr << "ERROR: errno=" << 1 << " opening file"
                << " at " << __FILE__ << ":" << __LINE__
                << std::endl << std::flush;
      exit(1);
  }

  unsigned nBlock = 20;
  double tmpa, tmpb, tmpc;
  vector<double > tt;
  vector<double > pp;

  while (fscanf (fp, "%lf%lf%lf", &tmpa, &tmpb, &tmpc) == 3){
    tt.push_back (tmpb);
    pp.push_back (tmpc);
  }

  
  BlockAverage ba;
  ba.processData (pp, nBlock);
  
  unsigned nDataInBlock = pp.size() / nBlock;
  BlockAverage_acc bacc (nDataInBlock);  
  for (unsigned ii = 0; ii < pp.size(); ++ii){
    bacc.deposite (pp[ii]);
  }
  bacc.calculate();

  cout << "ndata used: " << ba.getNumDataUsed() << endl;
  cout << "ndata used: " << bacc.getNumDataUsed() << endl;
  cout << "ba:     " << ba.getAvg() << " " << ba.getAvgError() << endl;
  cout << "ba_acc: " << bacc.getAvg() << " " << bacc.getAvgError() << endl;
  
  
}
