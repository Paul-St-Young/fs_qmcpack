//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file ECPComponentBuilderBuilder.h
 * @brief Declaration of a builder class for an ECP component for an ionic type
 */
#ifndef QMCPLUSPLUS_ECPCOMPONENT_BUILDER_H
#define QMCPLUSPLUS_ECPCOMPONENT_BUILDER_H
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/LocalECPotential.h"
#include "QMCHamiltonians/NonLocalECPotential.h"

namespace qmcplusplus
{

struct ECPComponentBuilder: public MPIObjectBase, public QMCTraits
{

  typedef LocalECPotential::GridType GridType;
  typedef ParticleSet::Scalar_t mRealType;
  typedef OneDimGridBase<mRealType> mGridType;
  typedef LocalECPotential::RadialPotentialType RadialPotentialType;

  int NumNonLocal;
  int Lmax, Llocal, Nrule;
  RealType Zeff;
  RealType RcutMax;
  std::string Species;
  mGridType *grid_global;
  std::map<std::string,mGridType*> grid_inp;
  RadialPotentialType* pp_loc;
  NonLocalECPComponent* pp_nonloc;
  std::map<std::string,int> angMon;

  ECPComponentBuilder(const std::string& aname, Communicate* c);

  bool parse(const std::string& fname, xmlNodePtr cur);
  bool put(xmlNodePtr cur);
  void addSemiLocal(xmlNodePtr cur);
  void buildLocal(xmlNodePtr cur);
  void buildSemiLocalAndLocal(std::vector<xmlNodePtr>& semiPtr);

  bool parseCasino(const std::string& fname, xmlNodePtr cur); //std::string& fname, RealType rc);
  //bool parseCasino(std::string& fname, RealType rc);
  // This sets the spherical quadrature rule used to apply the
  // projection operators.  rule can be 1 to 7.  See
  // J. Chem. Phys. 95 (3467) (1991)
  // Rule     # points     lexact
  //  1           1          0
  //  2           4          2
  //  3           6          3
  //  4          12          5
  //  5          18          5
  //  6          26          7
  //  7          50         11
  void SetQuadratureRule(int rule);
  void CheckQuadratureRule(int lexact);

  mGridType* createGrid(xmlNodePtr cur, bool useLinear=false);
  RadialPotentialType* createVrWithBasisGroup(xmlNodePtr cur, mGridType* agrid);
  RadialPotentialType* createVrWithData(xmlNodePtr cur, mGridType* agrid, int rCorrection=0);

  void doBreakUp(const std::vector<int>& angList, const Matrix<mRealType>& vnn,
                 RealType rmax, mRealType Vprefactor=1.0);

  void printECPTable();
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: yeluo $
 * $Revision: 7167 $   $Date: 2016-10-07 16:46:07 -0600 (Fri, 07 Oct 2016) $
 * $Id: ECPComponentBuilder.h 7167 2016-10-07 22:46:07Z yeluo $
 ***************************************************************************/
