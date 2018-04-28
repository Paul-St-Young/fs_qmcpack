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
    
    



#include "Utilities/OhmmsSpecies.h"

SpeciesBase::SpeciesBase()
{
  TotalNum = 0;
  Name.reserve(10); // expect less than 10 species
}

SpeciesBase::~SpeciesBase()
{
  AttribList_t::iterator dit = d_attrib.begin();
  for(; dit != d_attrib.end(); ++dit)
    delete (*dit);
}

int SpeciesBase::addAttrib()
{
  int n = d_attrib.size();
  d_attrib.push_back(new SpeciesAttrib_t(TotalNum));
  return n;
}

/***************************************************************************
 * $RCSfile$   $Author: abenali $
 * $Revision: 7138 $   $Date: 2016-09-27 17:45:29 -0600 (Tue, 27 Sep 2016) $
 * $Id: OhmmsSpecies.cpp 7138 2016-09-27 23:45:29Z abenali $
 ***************************************************************************/


/***************************************************************************
 * $RCSfile$   $Author: abenali $
 * $Revision: 7138 $   $Date: 2016-09-27 17:45:29 -0600 (Tue, 27 Sep 2016) $
 * $Id: OhmmsSpecies.cpp 7138 2016-09-27 23:45:29Z abenali $
 ***************************************************************************/

