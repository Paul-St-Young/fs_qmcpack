//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_DUMMY_H
#define QMCPLUSPLUS_DUMMY_H
#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers
 *@brief A dummy QMCDriver for testing
 */
class DummyQMC: public QMCDriver
{
public:
  /// Constructor.
  DummyQMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  bool run();
  bool put(xmlNodePtr cur);

private:
  /// Copy Constructor (disabled)
  DummyQMC(const DummyQMC& a): QMCDriver(a) { }
  /// Copy operator (disabled).
  DummyQMC& operator=(const DummyQMC&)
  {
    return *this;
  }

};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: abenali $
 * $Revision: 7138 $   $Date: 2016-09-27 17:45:29 -0600 (Tue, 27 Sep 2016) $
 * $Id: DummyQMC.h 7138 2016-09-27 23:45:29Z abenali $
 ***************************************************************************/
