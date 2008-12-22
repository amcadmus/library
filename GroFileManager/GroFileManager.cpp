#include "GroFileManager.h"
#include <iterator>

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
  
//   fprintf (fp, "\n%d\n", resdindex.size());
//   for (int i = 0; i < resdindex.size(); ++i){
//     fprintf (fp, "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
// 	     resdindex[i], (char *)(resdname[i].c_str()), (char *)(atomname[i].c_str()), atomindex[i], 
// 	     posi[i][0], posi[i][1], posi[i][2],
// 	     velo[i][0], velo[i][1], velo[i][2]);
//   }
//   fprintf (fp, "%f %f %f\n", boxsize[0], boxsize[1], boxsize[2]);
  fclose (fp);
}


int main ()
{
  std::string  name ;
  std::vector<int >  resdindex;
  std::vector<std::string >   resdname;
  std::vector<std::string >   atomname;
  std::vector<int >  atomindex;
  std::vector<std::vector<double > >  posi;
  std::vector<std::vector<double > >  velo;
  std::vector<double >  boxsize;
			   
  GroFileManager::read ("n32000_size40.0_seed67321.gro",
			resdindex, resdname, atomname, atomindex,
			posi, velo, boxsize);

  GroFileManager::write ("tmp.gro",
			 resdindex, resdname, atomname, atomindex,
			 posi, velo, boxsize);

  FILE * fp = fopen ("tmp", "w");
  fprintf (fp, "%d\n", resdindex.size());
  for (int i =0; i < int(resdindex.size()); ++i){
    fprintf (fp, "%8.3f%8.3f%8.3f\n", posi[i][0], posi[i][1], posi[i][2]);
  }
  fprintf (fp, "%f\n", boxsize[0]);
  fclose (fp);
  
  return 0;
}
