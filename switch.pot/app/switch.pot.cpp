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
#include "StringSplit.h"

#define MaxLineLength 2048

namespace po = boost::program_options;
using namespace std;

int main(int argc, char * argv[])
{
  std::string ifile, ofile;
  double x0, x1, ;

  po::options_description desc ("Allow options");
  desc.add_options()
      ("help,h", "print this message")
      ("start,b", po::value<double > (&x0)->default_value (0.9), "start switch")
      ("end,e",   po::value<double > (&x1)->default_value (1.1), "end switch")
      ("input,f",  po::value<string > (&ifile)->default_value ("table.xvg"),   "the input table potential")
      ("output,o", po::value<string > (&ofile)->default_value ("out.xvg"), "the output table potential");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify (vm);
  if (vm.count("help")){
    std::cout << desc<< "\n";
    return 0;
  }

  ifstream fpname (ifile.c_str());
  if (!fpname){
    std::cerr << "cannot open file " << ifile << std::endl;
    return 1;
  }
  while (fpname.getline(nameline, MaxLineLength)){

  
  return 0;
}
