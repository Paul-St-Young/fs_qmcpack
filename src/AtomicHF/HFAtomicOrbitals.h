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
    
    


#ifndef OHMMS_ATOMICHARTREEFOCK_TYPES_H
#define OHMMS_ATOMICHARTREEFOCK_TYPES_H

#include "AtomicHF/YlmRnlSet.h"

/**@namespace ohmmshf
 *@brief Define basic data types for the applications.
 * In order to reduce complier-time complexity and to enable switching
 * between  libraries for array and expression template,
 * basic data types are defined.
 */
namespace ohmmshf
{
typedef YlmRnlSet<double> HFAtomicOrbitals;
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: abenali $
 * $Revision: 7138 $   $Date: 2016-09-27 17:45:29 -0600 (Tue, 27 Sep 2016) $
 * $Id: HFAtomicOrbitals.h 7138 2016-09-27 23:45:29Z abenali $
 ***************************************************************************/
