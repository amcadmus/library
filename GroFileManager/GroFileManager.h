#ifndef __GroFileManager_wanghan__
#define __GroFileManager_wanghan__

#include <vector>
#include <iostream>
#include <string>

namespace GroFileManager{
    void read (const std::string & name ,
	       std::vector<int > & resdindex,
	       std::vector<std::string > &  resdname,
	       std::vector<std::string > &  atomname,
	       std::vector<int > & atomindex,
	       std::vector<std::vector<double > > & posi,
	       std::vector<std::vector<double > > & velo,
	       std::vector<double > & boxsize);
    void write (const std::string & name ,
		const std::vector<int > & resdindex,
		const std::vector<std::string > &  resdname,
		const std::vector<std::string > &  atomname,
		const std::vector<int > & atomindex,
		const std::vector<std::vector<double > > & posi,
		const std::vector<std::vector<double > > & velo,
		const std::vector<double > & boxsize);
};









#endif
