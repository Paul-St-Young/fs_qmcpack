//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_VECTORREF_H
#define OHMMS_VECTORREF_H

template<class T>
struct VectorRef
{

  typedef T value_type;
  VectorRef(T* datain):dptr(datain) {}

  inline T& operator[](int i)
  {
    return dptr[i];
  }
  inline T operator[](int i) const
  {
    return dptr[i];
  }
  T* dptr;
};

#endif

/***************************************************************************
 * $RCSfile$   $Author: abenali $
 * $Revision: 7138 $   $Date: 2016-09-27 17:45:29 -0600 (Tue, 27 Sep 2016) $
 * $Id: OhmmsVectorRef.h 7138 2016-09-27 23:45:29Z abenali $
 ***************************************************************************/

