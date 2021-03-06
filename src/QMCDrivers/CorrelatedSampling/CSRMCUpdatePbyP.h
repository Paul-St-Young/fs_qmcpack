//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_CS_VMC_UPDATEPBYP_H
#define QMCPLUSPLUS_CS_VMC_UPDATEPBYP_H
#include "QMCDrivers/CorrelatedSampling/CSUpdateBase.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers MultiplePsi ParticleByParticle
 * @brief Implements the VMC algorithm
 */
class CSVMCUpdatePbyP: public CSUpdateBase
{

public:
  /// Constructor.
  CSVMCUpdatePbyP(MCWalkerConfiguration& w,
                  TrialWaveFunction& psi,
                  QMCHamiltonian& h,
                  RandomGenerator_t& rg);

  ~CSVMCUpdatePbyP();

  void advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure);

private:
};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1592 $   $Date: 2007-01-04 16:48:00 -0600 (Thu, 04 Jan 2007) $
 * $Id: CSVMCUpdatePbyP.h 1592 2007-01-04 22:48:00Z jnkim $
 ***************************************************************************/
