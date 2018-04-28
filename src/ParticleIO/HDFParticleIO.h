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
    
    


#ifndef OHMMS_PARTICLE_INPUTOUTPUT_HDF_UTILITY_H
#define OHMMS_PARTICLE_INPUTOUTPUT_HDF_UTILITY_H

#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/RecordProperty.h"
#include "Particle/ParticleSet.h"
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus
{

class HDFParticleParser: public ParticleTags
{

public:

  typedef ParticleSet Particle_t;

  HDFParticleParser(Particle_t& aptcl):ref_(aptcl) { }

  ///reading from a file
  bool put(const char*);

  bool put(xmlNodePtr cur);

private:
  Particle_t& ref_;
};

class HDFSaveParticle:
  public ParticleTags,
  public RecordProperty
{

public:

  typedef ParticleSet Particle_t;

  HDFSaveParticle(Particle_t& pin): ref_(pin) { }

  ~HDFSaveParticle();

  void reset(const char* fileroot, bool append=false);

  void report(int iter);

  void finalize() { }

  bool put(xmlNodePtr cur);

private:

  Particle_t& ref_;
  std::string FileRoot;

};
}

#endif

/***************************************************************************
 * $RCSfile$   $Author: abenali $
 * $Revision: 7138 $   $Date: 2016-09-27 17:45:29 -0600 (Tue, 27 Sep 2016) $
 * $Id: HDFParticleIO.h 7138 2016-09-27 23:45:29Z abenali $
 ***************************************************************************/
