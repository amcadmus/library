#include "GroFileManager.h"
#include <iterator>
#include <iostream>
#include <fstream>

bool GroFileManager::readTop (const std::string & filename,
			      std::vector<std::string > & molnames,
			      std::vector<int > & nmols)
{
  molnames.clear();
  nmols.clear();

  std::ifstream in(filename.c_str());
  if (in.bad()){
    std::cerr << "cannot open file " << filename << std::endl;
    return false;
  }
  char line [1024];
  std::string target ("[ molecules ]");
  bool find = false;
  while (!in.eof()){
    in.getline (line, 1024, '\n');
    if (target == std::string (line)) {
      find = true;
      break;
    }
  }
  if (!find){
    std::cerr << "cannot find [ molecules ] in file " << filename 
	      << ". please check there is no space after \"]\"\n";
    return false;
  }
  
//   while (!(in.getline (line, 1024, '\n')).eof()){
// //     if (line[0] == '['){
// //       break;
// //     }
// //     char name[1024];
// //     int number;
// //     sscanf (line, "%s%d", name, &number);
// //     molnames.push_back (std::string(name));
// //     nmols.push_back (number);
//   }
  

  std::string name;
  int number;
  while ( !(in >> name >> number).eof()){
    if (name[0] == '['){
      break;
    }
    if (name.empty()){
      break;
    }
//     std::cout << name << std::endl;
    molnames.push_back (name);
    nmols.push_back (number);
  }
  

  return true;
}



template <typename UnitaryFunction1, typename UnitaryFunction2,
	  typename UnitaryFunction3, typename UnitaryFunction4,
	  typename UnitaryFunction5, typename UnitaryFunction6>
bool GroFileManager::writePotenFile (const double & rmin, const double & rcut, 
				     const double & interval,
				     UnitaryFunction1 & f, UnitaryFunction2 & fp,
				     UnitaryFunction3 & g, UnitaryFunction4 & gp,
				     UnitaryFunction5 & h, UnitaryFunction6 & hp,
				     const std::string & filename)
{
  FILE * filep = fopen (filename.c_str(), "w");
  if(filep == NULL){
    std::cerr << "cannot open file " << filename << std::endl;
    return false;
  }

  double upper = rcut + 1;
  double nx;
  if ( int(upper / interval) != upper / interval)
    nx = int(upper / interval) + 1;
  else 
    nx = int(upper / interval);
  upper = interval * nx;
  
  int i = 0;
  for (i = 0; i <= nx + 1; ++i){
    double x = i * interval;
    if (x < rmin){
      fprintf (filep, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
	       x, 0., 0., 0., 0., 0., 0.);
    }
    else {
      fprintf (filep, "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
	       x, f(x), -fp(x), g(x), -gp(x), h(x), -hp(x));
    }
  }
  
  fclose (filep);
  return true;
}



void GroFileManager::read (const std::string & name ,
			   std::vector<int > & resdindex,
			   std::vector<std::string > &  resdname,
			   std::vector<std::string > &  atomname,
			   std::vector<int > & atomindex,
			   std::vector<std::vector<double > > & posi,
			   std::vector<std::vector<double > > & velo,
			   std::vector<double > & boxsize)
{
  FILE * fp = fopen (name.c_str(), "r");
  if (fp == NULL){
    std::cerr << "cannot open file " << name << std::endl;
    return;
  }

  resdindex.clear();
  resdname.clear();
  atomname.clear();
  atomindex.clear();
  posi.clear();
  velo.clear();
  boxsize.resize(3);
  
  while (fgetc(fp) != '\n');
  int npart;
  int n1, n2;
  char s1[10], s2[10];
  double a, b, c;
  double d, e, f;
  fscanf (fp, "%d", &npart);
  for (int i = 0; i < npart; ++ i){
    fscanf (fp, "%5d%5s%5s%5d%lf%lf%lf%lf%lf%lf",
	    &n1, s1, s2, &n2, &a, &b, &c, &d, &e, &f);
    resdindex.push_back (n1);
    resdname.push_back (std::string(s1));
    atomname.push_back (std::string(s2));
    atomindex.push_back (n2);
    std::vector<double > tmp;
    tmp.push_back(a);
    tmp.push_back(b);
    tmp.push_back(c);
    posi.push_back(tmp);
    tmp.clear();
    tmp.push_back(d);
    tmp.push_back(e);
    tmp.push_back(f);
    velo.push_back(tmp);    
  }
  fscanf (fp, "%8lf%8lf%8lf", &boxsize[0], &boxsize[1], &boxsize[2]);
  
  fclose (fp);

}

void GroFileManager::write (const std::string & name ,
			    const std::vector<int > & resdindex,
			    const std::vector<std::string > &  resdname,
			    const std::vector<std::string > &  atomname,
			    const std::vector<int > & atomindex,
			    const std::vector<std::vector<double > > & posi,
			    const std::vector<std::vector<double > > & velo,
			    const std::vector<double > & boxsize)
{
  FILE * fp = fopen(name.c_str(), "w");
  if (fp == NULL){
    std::cerr << "cannot open file " << name << std::endl;
    return;
  }
  std::copy (atomname.begin(), atomname.end(), std::ostream_iterator<std::string>(std::cout, "\n"));
  
  fprintf (fp, "\n%d\n", resdindex.size());
  for (int i = 0; i < int(resdindex.size()); ++i){
    fprintf (fp, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
	     resdindex[i], (char *)(resdname[i].c_str()), (char *)(atomname[i].c_str()), atomindex[i], 
	     posi[i][0], posi[i][1], posi[i][2],
	     velo[i][0], velo[i][1], velo[i][2]);
  }
  fprintf (fp, "%f %f %f\n", boxsize[0], boxsize[1], boxsize[2]);

  fclose (fp);
}


struct F 
{
  double operator () (double x)
      {
	return 1./x;
      }
}
    ;
struct Zero
{
  double operator () (double x)
      {
	return 0;
      }
}
    ;

