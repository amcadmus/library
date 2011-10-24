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
#include <boost/program_options.hpp>

const double kB = 1.38065e-23; // J / K
const double avogadro = 6.0221415e23; // mol-1

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  double c6, c12;
  double temperature;

  cout.precision (4);
  cout << (scientific);
  
  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("c6",  po::value<double > (&c6) ->default_value(2.6171e-3),  "c6  [kJ/(mol nm^6)]")
      ("c12", po::value<double > (&c12)->default_value(2.6331e-6),  "c12 [kJ/(mol nm^12)]")
      ("temperature,t", po::value<double > (&temperature)->default_value(300),  "temperature [K]");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  printf ("# c6 \t\t: %.6e [kJ/(mol nm^6)]\n", c6);
  printf ("# c12\t\t: %.6e [kJ/(mol nm^12)]\n", c12);
  printf ("# T  \t\t: %.2f [K]\n", temperature);

  double sigma = pow (c12 / c6, 1./6.);
  double epsilon = c6 / (4. * sigma * sigma * sigma * sigma * sigma * sigma);
  double kT = kB * avogadro * temperature * 1e-3;
  
  printf ("epsilon\t\t: %.6f [kJ/mol]\n", epsilon);
  printf ("sigma  \t\t: %.6f [nm]\n", sigma);
  printf ("kT     \t\t: %.6f [kJ/mol]\n", kT);
  printf ("epsilon / kT\t: %.6f \n", epsilon / kT);
  

  return 0;
}
