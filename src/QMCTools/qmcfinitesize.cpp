//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "QMCApp/ParticleSetPool.h"

#include "QMCApp/QMCAppBase.h"
#include "QMCTools/QMCFiniteSize/QMCFiniteSize.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "QMCTools/QMCFiniteSize/SkParserASCII.h"
#include "QMCTools/QMCFiniteSize/SkParserScalarDat.h"

#include "Numerics/OneDimGridBase.h"

//Purpose of this routine is to compute the finite size effects
//for a given simulation cell from post-processed QMC Data.  For
//the potential, this is done by splining the structure factor
//and performing the integral given in Holzmann et al., PRB 035126 (2016)
//
//Input:  Cell geometry.  Ion positions, cell geometries, etc, are taken from main.xml file.
//                        Code recognizes ESHDF5 declarations in main.xml.
//        S(k):  This is the electron-electron fluctuation structure factor.
//
//Returns: (E(N=infty)-E(N)) for the given simulation cell.

using namespace qmcplusplus;
typedef QMCTraits::RealType RealType;
typedef QMCTraits::PosType  PosType;
typedef SkParserBase::Grid_t Grid_t;

int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc,argv);
  OhmmsInfo Welcome("qmcfinitesize",OHMMS::Controller->rank());
  Random.init(0,1,-1);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.setf(std::ios::right,std::ios::adjustfield);
  std::cout.precision(12);
 
  std::cout<<"Usage:  qmcfinitesize [main.xml] --[skformat] [SK_FILE]\n"; 
  std::cout<<"  --ascii:  assume S(k) is given in kx ky kz sk sk_err format.  Header necessary.\n";
  std::cout<<"  --scalardat:  parses 'skall' observable from files with energy.pl output format.\n";
  
  SkParserBase* skparser(NULL);
  int iargc=2;
  
  while(iargc+1<argc)
  {
    std::string a(argv[iargc]);
    std::string anxt(argv[iargc+1]);
    std::cout<<" "<<a<<"  "<<anxt<<std::endl;
    if(a=="--ascii")
    {
      skparser=new SkParserASCII();
      skparser->parse(anxt);
    }
    else if(a=="--scalardat")
    {
      skparser=new SkParserScalarDat();
      skparser->parse(anxt);
    }
    iargc++;
  }
  
  if(skparser==NULL) APP_ABORT("qmcfinitesize:  skparser failed to initialize");
  
  //We are only going to do one main.xml file for now.  
  std::cout<<"Checking to see if "<<argv[1]<<" is a valid xml file\n";
  std::string a(argv[1]);
  std::size_t found = a.find(".xml");
  if(found==std::string::npos) APP_ABORT("qmcfinitesize: qmcinput file not an .xml file.");
  std::cout<<"Initializing QMCFiniteSize with skparser\n";
  QMCFiniteSize qmcfs(skparser);
  std::cout<<"Calling parse with "<<argv[1]<<" \n";
  qmcfs.parse(std::string(argv[1]));
  qmcfs.validateXML();
  qmcfs.execute(); 
  //LinearGrid<RealType> grid_x;
  //LinearGrid<RealType> grid_y;
  //LinearGrid<RealType> grid_z;
 // Grid_t grid_x;
 // Grid_t grid_y;
 // Grid_t grid_z;
 // cout<<"Calling skparser->get_grid(x,y,z)\n"; 
 // skparser->get_grid(grid_x,grid_y,grid_z);
 // cout<<"Done\n";
  
 
 // QMCFiniteSize chiesa;
 // chiesa.execute(); 
  OHMMS::Controller->finalize();
  delete skparser;
  
  return 0;
}

