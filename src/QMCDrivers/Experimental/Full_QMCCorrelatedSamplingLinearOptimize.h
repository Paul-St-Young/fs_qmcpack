//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file QMCCSLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCCSLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCCSLINEAROPTIMIZATION_VMCSINGLE_H

#include "QMCDrivers/QMCLinearOptimize.h"
#include "QMCDrivers/VMC/VMCLinearOptOMP.h"
#include "QMCDrivers/QMCCSLinearOptimizeWFmanagerOMP.h"
#include "Optimize/OptimizeBase.h"

namespace qmcplusplus
{

///forward declaration of a cost function
class QMCCostFunctionBase;
class HamiltonianPool;

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCCSLinearOptimize: public QMCLinearOptimize
{
public:

  ///Constructor.
  QMCCSLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                      QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool);

  ///Destructor
  ~QMCCSLinearOptimize();

  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);

private:
  VMCLinearOptOMP* vmcCSEngine;
  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  /// switch to control whether NRCOptimization::lineoptimization() is used or somethign else
  std::string MinMethod, GEVtype;

  RealType stabilizerScale, bigChange, exp0, stepsize;
  RealType Lambda;
  int nstabilizers;
  /// percent variance or H2 to mix in
  RealType w_beta;

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: QMCCSLinearOptimize.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
