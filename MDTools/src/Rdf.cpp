#include "Rdf.h"
#include <iostream>
#include <stdio.h>

void Rdf::
reinit (const ValueType rup_,
	const ValueType refh_,
	const ValueType x0_,
	const ValueType x1_,
	const int nDataInBlock)
{
  rup = rup_;
  binSize = refh_;
  x0 = x0_;
  x1 = x1_;
  if (x0 > x1) {
    ValueType tmpx = x0;
    x0 = x1;
    x1 = tmpx;
  }  
  nbins = rup / refh_ + 1;
  rup = binSize * (nbins - .5);
  offset = .5 * binSize;
  hist.resize (nbins);
  value.resize (nbins);
  error.resize (nbins);
  std::fill (hist.begin(), hist.end(), 0.);
  std::fill (value.begin(), value.end(), 0.);
  std::fill (error.begin(), error.end(), 0.);
  avg.resize (nbins);
  for (unsigned ii = 0; ii < avg.size(); ++ii){
    avg[ii].reinit (nDataInBlock);
  }
  nframe = 0;
  rho = 0.;
  natom = 0.;
}

// void Rdf::
// deposit (const std::vector<std::vector<ValueType> > & coord,
// 	 const VectorType & box,
// 	 const CellList & clist)
// {
//   int xiter = rup / clist.getCellSize().x;
//   if (xiter * clist.getCellSize().x < rup) xiter ++;
//   int yiter = rup / clist.getCellSize().y;
//   if (yiter * clist.getCellSize().y < rup) yiter ++;
//   int ziter = rup / clist.getCellSize().z;
//   if (ziter * clist.getCellSize().z < rup) ziter ++;

//   IntVectorType nCell = clist.getNumCell();
//   ValueType myNatom = 0.;
  
//   for (int ix = 0; ix < nCell.x; ++ix){
//     for (int iy = 0; iy < nCell.y; ++iy){
//       for (int iz = 0; iz < nCell.z; ++iz){
// 	unsigned iid = clist.index3to1 (ix, iy, iz);
// 	for (unsigned ii = 0; ii < clist.getList()[iid].size(); ++ii){
// 	  VectorType icoord;
// 	  icoord.x = coord[clist.getList()[iid][ii]][0];
// 	  if (x1 != 0. && (!(icoord.x >= x0 && icoord.x < x1))) continue;
// 	  icoord.y = coord[clist.getList()[iid][ii]][1];
// 	  icoord.z = coord[clist.getList()[iid][ii]][2];
// 	  myNatom += 1.;
// 	  for (int dx = -xiter; dx <= xiter; ++dx){
// 	    int jx = ix + dx;
// 	    if (jx < 0) jx += nCell.x;
// 	    else if (jx >= nCell.x) jx -= nCell.x;
// 	    for (int dy = -yiter; dy <= yiter; ++dy){
// 	      int jy = iy + dy;
// 	      if (jy < 0) jy += nCell.y;
// 	      else if (jy >= nCell.y) jy -= nCell.y;
// 	      for (int dz = -ziter; dz <= ziter; ++dz){
// 		int jz = iz + dz;
// 		if (jz < 0) jz += nCell.z;
// 		else if (jz >= nCell.z) jz -= nCell.z;
// 		unsigned jid = clist.index3to1 (jx, jy, jz);
// 		bool sameCell (0 == dx && 0 == dy && 0 == dz);
// 		for (unsigned jj = 0; jj < clist.getList()[jid].size(); ++jj){
// 		  if (sameCell && ii == jj) continue;
// 		  VectorType jcoord;
// 		  jcoord.x = coord[clist.getList()[jid][jj]][0];
// 		  jcoord.y = coord[clist.getList()[jid][jj]][1];
// 		  jcoord.z = coord[clist.getList()[jid][jj]][2];
// 		  VectorType diff;
// 		  diff.x = icoord.x - jcoord.x;
// 		  diff.y = icoord.y - jcoord.y;
// 		  diff.z = icoord.z - jcoord.z;
// 		  if      (diff.x < -.5 * box.x) diff.x += box.x;
// 		  else if (diff.x >= .5 * box.x) diff.x -= box.x;
// 		  if      (diff.y < -.5 * box.y) diff.y += box.y;
// 		  else if (diff.y >= .5 * box.y) diff.y -= box.y;
// 		  if      (diff.z < -.5 * box.z) diff.z += box.z;
// 		  else if (diff.z >= .5 * box.z) diff.z -= box.z;
// 		  ValueType dr = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
// 		  dr = sqrt (dr);
// 		  unsigned index = (dr + offset) / binSize;
// 		  if (dr < rup){
// 		    if (index >= unsigned(nbins)){
// 		      // printf ("# dr: %f, index: %d, rup: %f, nbins: %d\n",
// 		      // 	      dr, index, rup, nbins);
// 		      index = nbins - 1;
// 		    }
// 		    hist[index] += 1.;
// 		  }
// 		}
// 	      }	    
// 	    }
// 	  }
// 	}
//       }
//     }
//   }

//   nframe ++;
//   if (x1 == x0){
//     rho += myNatom / (box.x * box.y * box.z);
//   }
//   else {
//     rho += myNatom / ((x1 - x0) * box.y * box.z);
//   }
//   natom += myNatom;
// }



void Rdf::
deposit (const std::vector<std::vector<ValueType> > & coord,
	 const VectorType & box,
	 const CellList & clist)
{
  int xiter = rup / clist.getCellSize().x;
  if (xiter * clist.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist.getCellSize().y;
  if (yiter * clist.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist.getCellSize().z;
  if (ziter * clist.getCellSize().z < rup) ziter ++;

  IntVectorType nCell = clist.getNumCell();
  ValueType myNatom = 0.;
  std::vector<ValueType > tmphist (hist.size(), 0.);
  
  for (unsigned iCellIndex = 0;
       iCellIndex < unsigned(nCell.x * nCell.y * nCell.z);
       ++iCellIndex){
    std::vector<unsigned > neighborCellIndex =
	clist.neighboringCellIndex (iCellIndex, IntVectorType (xiter, yiter, ziter));
    for (unsigned iNeighborCellIndex = 0;
	 iNeighborCellIndex < neighborCellIndex.size();
	 ++iNeighborCellIndex){
      unsigned jCellIndex = neighborCellIndex[iNeighborCellIndex];
      for (unsigned ii = 0; ii < clist.getList()[iCellIndex].size(); ++ii){
	VectorType icoord;
	icoord.x = coord[clist.getList()[iCellIndex][ii]][0];
	if (x1 != 0. && (!(icoord.x >= x0 && icoord.x < x1))) continue;
	if (iNeighborCellIndex == 0) myNatom += 1.;
	icoord.y = coord[clist.getList()[iCellIndex][ii]][1];
	icoord.z = coord[clist.getList()[iCellIndex][ii]][2];
	bool sameCell (iCellIndex == jCellIndex);
	for (unsigned jj = 0; jj < clist.getList()[jCellIndex].size(); ++jj){
	  if (sameCell && ii == jj) continue;	    
	  VectorType jcoord;
	  jcoord.x = coord[clist.getList()[jCellIndex][jj]][0];
	  jcoord.y = coord[clist.getList()[jCellIndex][jj]][1];
	  jcoord.z = coord[clist.getList()[jCellIndex][jj]][2];
	  VectorType diff;
	  diff.x = - icoord.x + jcoord.x;
	  diff.y = - icoord.y + jcoord.y;
	  diff.z = - icoord.z + jcoord.z;
	  if      (diff.x < -.5 * box.x) diff.x += box.x;
	  else if (diff.x >= .5 * box.x) diff.x -= box.x;
	  if      (diff.y < -.5 * box.y) diff.y += box.y;
	  else if (diff.y >= .5 * box.y) diff.y -= box.y;
	  if      (diff.z < -.5 * box.z) diff.z += box.z;
	  else if (diff.z >= .5 * box.z) diff.z -= box.z;
	  ValueType dr = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
	  dr = sqrt (dr);
	  unsigned index = (dr + offset) / binSize;
	  if (dr < rup){
	    if (index >= unsigned(nbins)){
	      // printf ("# dr: %f, index: %d, rup: %f, nbins: %d\n",
	      // 	      dr, index, rup, nbins);
	      index = nbins - 1;
	    }
	    tmphist[index] += 1.;
	  }
	}
      }
    }
  }

  // printf ("\n");
  nframe ++;
  if (x1 == x0){
    rho += myNatom / (box.x * box.y * box.z);
  }
  else {
    rho += myNatom / ((x1 - x0) * box.y * box.z);
  }
  natom += myNatom;

  for (unsigned ii = 0; ii < avg.size(); ++ii){
    avg[ii].deposite (tmphist[ii]/double(myNatom));
    hist[ii] += tmphist[ii];
  }
}



void Rdf::
deposit (const std::vector<std::vector<ValueType> > & coord1,
	 const std::vector<std::vector<ValueType> > & coord2,
	 const VectorType & box,
	 const CellList & clist1,
	 const CellList & clist2)
{
  int xiter = rup / clist2.getCellSize().x;
  if (xiter * clist2.getCellSize().x < rup) xiter ++;
  int yiter = rup / clist2.getCellSize().y;
  if (yiter * clist2.getCellSize().y < rup) yiter ++;
  int ziter = rup / clist2.getCellSize().z;
  if (ziter * clist2.getCellSize().z < rup) ziter ++;

  IntVectorType nCell = clist1.getNumCell();
  ValueType myNatom1 = 0.;
  ValueType myNatom2 = 0.;
  for (unsigned ii = 0; ii < clist2.getTotalNumCell(); ++ii){
    myNatom2 += clist2.getList()[ii].size();
  }
  std::vector<ValueType > tmphist (hist.size(), 0.);
  
  for (unsigned iCellIndex = 0;
       iCellIndex < unsigned(nCell.x * nCell.y * nCell.z);
       ++iCellIndex){
    std::vector<unsigned > neighborCellIndex =
	clist2.neighboringCellIndex (iCellIndex, IntVectorType (xiter, yiter, ziter));
    for (unsigned iNeighborCellIndex = 0;
	 iNeighborCellIndex < neighborCellIndex.size();
	 ++iNeighborCellIndex){
      unsigned jCellIndex = neighborCellIndex[iNeighborCellIndex];
      for (unsigned ii = 0; ii < clist1.getList()[iCellIndex].size(); ++ii){
	VectorType icoord;
	icoord.x = coord1[clist1.getList()[iCellIndex][ii]][0];
	if (x1 != 0. && (!(icoord.x >= x0 && icoord.x < x1))) continue;
	if (iNeighborCellIndex == 0) myNatom1 += 1.;
	icoord.y = coord1[clist1.getList()[iCellIndex][ii]][1];
	icoord.z = coord1[clist1.getList()[iCellIndex][ii]][2];
	// bool sameCell (iCellIndex == jCellIndex);
	for (unsigned jj = 0; jj < clist2.getList()[jCellIndex].size(); ++jj){
	  // if (sameCell && ii == jj) continue;	    
	  VectorType jcoord;
	  jcoord.x = coord2[clist2.getList()[jCellIndex][jj]][0];
	  jcoord.y = coord2[clist2.getList()[jCellIndex][jj]][1];
	  jcoord.z = coord2[clist2.getList()[jCellIndex][jj]][2];
	  VectorType diff;
	  diff.x = - icoord.x + jcoord.x;
	  diff.y = - icoord.y + jcoord.y;
	  diff.z = - icoord.z + jcoord.z;
	  if      (diff.x < -.5 * box.x) diff.x += box.x;
	  else if (diff.x >= .5 * box.x) diff.x -= box.x;
	  if      (diff.y < -.5 * box.y) diff.y += box.y;
	  else if (diff.y >= .5 * box.y) diff.y -= box.y;
	  if      (diff.z < -.5 * box.z) diff.z += box.z;
	  else if (diff.z >= .5 * box.z) diff.z -= box.z;
	  ValueType dr = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
	  dr = sqrt (dr);
	  unsigned index = (dr + offset) / binSize;
	  if (dr < rup){
	    if (index >= unsigned(nbins)){
	      // printf ("# dr: %f, index: %d, rup: %f, nbins: %d\n",
	      // 	      dr, index, rup, nbins);
	      index = nbins - 1;
	    }
	    tmphist[index] += 1.;
	  }
 	}
      }
    }
  }

  // printf ("\n");
  nframe ++;
  // if (x1 == x0){
  //   rho += myNatom2 / (box.x * box.y * box.z);
  // }
  // else {
  //   rho += myNatom2 / ((x1 - x0) * box.y * box.z);
  // }
  rho += myNatom2 / (box.x * box.y * box.z);
  natom += myNatom1;
  for (unsigned ii = 0; ii < avg.size(); ++ii){
    avg[ii].deposite (tmphist[ii]/double(myNatom1));
    hist[ii] += tmphist[ii];
  }
}



void Rdf::
calculate()
{
  rho /= double(nframe);
  for (int i = 0; i < nbins; ++i){
    hist[i] /= double(natom);
    avg[i].calculate();
  }
  {
    double r = 0.5 * binSize;
    hist[0] /= 4. / 3. * M_PI * r * r * r * rho;
    value[0] = avg[0].getAvg() / (4. / 3. * M_PI * r * r * r * rho);
    error[0] = avg[0].getAvgError() / (4. / 3. * M_PI * r * r * r * rho);
  }
  for (int i = 1; i < nbins; ++i){
    double r0 = (i-0.5) * binSize;
    double r1 = (i+0.5) * binSize;
    // double r01 = i * binSize;
    hist[i] /= 4. / 3. * M_PI * (r1*r1*r1 - r0*r0*r0) * rho;
    value[i] = avg[i].getAvg() / (4. / 3. * M_PI * (r1*r1*r1 - r0*r0*r0) * rho);
    error[i] = avg[i].getAvgError() / (4. / 3. * M_PI * (r1*r1*r1 - r0*r0*r0) * rho);
    // hist[i] /= 4. * M_PI * r0 * r1 * (r1 - r0) * rho;
  }
}



