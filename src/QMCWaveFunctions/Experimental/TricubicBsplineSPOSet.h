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
    
    
/** @file TricubicBsplineSPOSet.h
 * @brief Define TricubicBsplineSPOSet<T,bool ORTHO, bool TRUNC>
 *
 * Use specialization
 * TricubicBsplineSPOSet<T,true,false> : orthorhombic unit cell
 * TricubicBsplineSPOSet<T,false,false> : non-orthorhombic unit cell
 * TricubicBsplineSPOSet<T,true,true> : orthorhombic unit cell with localized orbitals
 * TricubicBsplineSPOSet<T,false,true> : non-orthorhombic unit cell with localized orbitals
 */
#ifndef TRICUBIC_BSPLINE_SINGLEORBITALSET_WITHSPECIALIZATION_H
#define TRICUBIC_BSPLINE_SINGLEORBITALSET_WITHSPECIALIZATION_H

#include "Numerics/TricubicBsplineGrid.h"
#include "Optimize/VarList.h"
#include <map>

namespace qmcplusplus
{

template<typename T, bool ORTHO, bool TRUNC>
class TricubicBsplineSPOSet: public TricubicBsplineTraits<T> {};

/** TricubicBsplineSPOSet specialized for Orthorhombic cell and no truncation*/
template<typename T>
class TricubicBsplineSPOSet<T,true,false>: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;
  typedef typename std::map<int,const StorageType*>::iterator  IteratorType;

  /** default constructure
   */
  TricubicBsplineSPOSet()
  { }
  ~TricubicBsplineSPOSet()
  {
  }

  inline void setTwistAngle(const PosType& tangle)
  {
  }

  inline void setGrid(const GridType& knots)
  {
    bKnots=knots;
  }

  ///empty reset
  void resetParameters(VarRegistry<real_type>& vlist)
  {
    ///DO NOTHING FOR NOW
  }

  inline void setGrid(real_type xi, real_type xf,
                      real_type yi, real_type yf, real_type zi, real_type zf,
                      int nx, int ny, int nz,
                      bool interp=true, bool periodic=true,bool openend=true)
  {
    bKnots.setGrid(xi,xf,yi,yf,zi,zf,nx,ny,nz,interp,periodic,openend);
  }

  /** add a orbital
   * @param i index of the orbital
   * @param data input data
   * @param curP interpolated data
   */
  void add(int i, const PosType& c, const StorageType& data, StorageType* curP)
  {
    IteratorType pit(P.find(i));
    if(pit == P.end())
    {
      bKnots.Init(data,*curP);
      P[i]=curP;
    }
  }

  void add(int i,const PosType& c,  StorageType* curP)
  {
    IteratorType pit(P.find(i));
    if(pit == P.end())
    {
      P[i]=curP;
    }
  }

  template<typename PV>
  inline void evaluate(const PosType& r, PV& vals)
  {
    bKnots.Find(r[0],r[1],r[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      vals[(*pit).first]=bKnots.evaluate(*((*pit).second));
      ++pit;
    }
  }

  template<typename PV, typename GV>
  inline void
  evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      int j((*pit).first);
      vals[j]=bKnots.evaluate(*((*pit).second),grads[j],laps[j]);
      ++pit;
    }
  }

  template<typename PM, typename GM>
  inline void
  evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      int j((*pit).first);
      vals(j,i)=bKnots.evaluate(*((*pit).second),grads(i,j),laps(i,j));
      ++pit;
    }
  }

private:
  GridType bKnots;
  std::map<int,const StorageType*> P;
};

/** TricubicBsplineSPOSet specialized for non-Orthorhombic cell and no truncation
 */
template<typename T>
class TricubicBsplineSPOSet<T,false,false>: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;
  typedef typename std::map<int,const StorageType*>::iterator  IteratorType;

  //going to use GGt=dot(Lattice.G,transpose(Lattice.G))
  using TricubicBsplineTraits<T>::GGt;
  //going to use Lattice
  using TricubicBsplineTraits<T>::Lattice;

  /** default constructure
  */
  TricubicBsplineSPOSet() { }

  ~TricubicBsplineSPOSet() { }

  inline void setTwistAngle(const PosType& tangle)
  {
  }

  inline void setGrid(const GridType& knots)
  {
    bKnots=knots;
  }

  ///empty reset
  void resetParameters(VarRegistry<real_type>& vlist)
  {
    ///DO NOTHING FOR NOW
  }

  inline void setGrid(real_type xi, real_type xf,
                      real_type yi, real_type yf, real_type zi, real_type zf,
                      int nx, int ny, int nz,
                      bool interp=true, bool periodic=true,bool openend=true)
  {
    bKnots.setGrid(0.0,1.0,0.0,1.0,0.0,1.0,nx,ny,nz,interp,periodic,openend);
  }

  /** add a orbital
   * @param i index of the orbital
   * @param data input data
   * @param curP interpolated data
   */
  void add(int i, const PosType& c, const StorageType& data, StorageType* curP)
  {
    IteratorType pit(P.find(i));
    if(pit == P.end())
    {
      bKnots.Init(data,*curP);
      P[i]=curP;
    }
  }

  void add(int i,const PosType& c,  StorageType* curP)
  {
    IteratorType pit(P.find(i));
    if(pit == P.end())
    {
      P[i]=curP;
    }
  }

  template<typename PV>
  inline void evaluate(const PosType& r, PV& vals)
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.Find(ru[0],ru[1],ru[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      vals[(*pit).first]=bKnots.evaluate(*((*pit).second));
      ++pit;
    }
  }

  template<typename PV, typename GV>
  inline void
  evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
  {
    PosType ru(Lattice.toUnit(r));
    TinyVector<value_type,3> gu;
    Tensor<value_type,3> hess;
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      int j((*pit).first);
      vals[j]=bKnots.evaluate(*((*pit).second),gu,hess);
      grads[j]=dot(Lattice.G,gu);
      laps[j]=trace(hess,GGt);
      ++pit;
    }
  }

  template<typename PM, typename GM>
  inline void
  evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
  {
    PosType ru(Lattice.toUnit(r));
    TinyVector<value_type,3> gu;
    Tensor<value_type,3> hess;
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    IteratorType pit(P.begin()), pit_end(P.end());
    while(pit != pit_end)
    {
      int j((*pit).first);
      vals(j,i)=bKnots.evaluate(*((*pit).second),gu,hess);
      grads(i,j)=dot(Lattice.G,gu);
      laps(i,j)=trace(hess,GGt);
      ++pit;
    }
  }

private:
  GridType bKnots;
  std::map<int,const StorageType*> P;
};

/** TricubicBsplineSPOSet specialized for Orthorhombic cell and no truncation*/
template<typename T>
class TricubicBsplineSPOSet<T,true,true>: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;
  typedef typename std::map<int,const StorageType*>::iterator  IteratorType;

  using TricubicBsplineTraits<T>::Rcut2;
  std::vector<PosType> Centers;

  /** default constructure
   *
   * Set Rcut2 to a large number so that everything counts
   */
  TricubicBsplineSPOSet()
  {
    Rcut2=1e6;
  }

  ~TricubicBsplineSPOSet()
  { }

  inline void setTwistAngle(const PosType& tangle)
  { }

  inline void setGrid(const GridType& knots)
  {
    bKnots=knots;
  }

  ///empty reset
  void resetParameters(VarRegistry<real_type>& vlist)
  {
    ///DO NOTHING FOR NOW
  }

  inline void setGrid(real_type xi, real_type xf,
                      real_type yi, real_type yf, real_type zi, real_type zf,
                      int nx, int ny, int nz,
                      bool interp=true, bool periodic=true,bool openend=true)
  {
    bKnots.setGrid(xi,xf,yi,yf,zi,zf,nx,ny,nz,interp,periodic,openend);
  }

  /** add a orbital
   * @param i index of the orbital
   * @param data input data
   * @param curP interpolated data
   */
  void add(int i, const PosType& c, const StorageType& data, StorageType* curP)
  {
    bKnots.Init(data,*curP);
    Centers.push_back(c);
    P.push_back(curP);
  }

  void add(int i, const PosType& c, StorageType* curP)
  {
    //Already exists
    if(i<Centers.size())
      return;
    Centers.push_back(c);
    P.push_back(curP);
  }

  template<typename PV>
  inline void evaluate(const PosType& r, PV& vals)
  {
    bKnots.Find(r[0],r[1],r[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
        vals[j]=0.0;//numeric_limits<T>::epsilon();
      else
        vals[j]=bKnots.evaluate(*P[j]);
    }
  }

  template<typename PV, typename GV>
  inline void
  evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        vals[j]=0.0;//numeric_limits<T>::epsilon();
        grads[j]=0.0;
        laps[j]=0.0;
      }
      else
        vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
    }
  }

  template<typename PM, typename GM>
  inline void
  evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
  {
    bKnots.FindAll(r[0],r[1],r[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        vals(j,i)=0.0; //numeric_limits<T>::epsilon();
        grads(i,j)=0.0;
        laps(i,j)=0.0;
      }
      else
        vals(j,i)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
    }
  }

private:
  GridType bKnots;
  std::vector<const StorageType*> P;
};

/** TricubicBsplineSPOSet specialized for non-Orthorhombic cell and truncation*/
template<typename T>
class TricubicBsplineSPOSet<T,false,true>: public TricubicBsplineTraits<T>
{
public:
  typedef typename TricubicBsplineTraits<T>::real_type   real_type;
  typedef typename TricubicBsplineTraits<T>::value_type  value_type;
  typedef typename TricubicBsplineTraits<T>::PosType     PosType;
  typedef typename TricubicBsplineTraits<T>::GridType    GridType;
  typedef typename TricubicBsplineTraits<T>::StorageType StorageType;
  typedef typename std::map<int,const StorageType*>::iterator  IteratorType;

  using TricubicBsplineTraits<T>::Rcut2;
  //going to use GGt=dot(Lattice.G,transpose(Lattice.G))
  using TricubicBsplineTraits<T>::GGt;
  //going to use Lattice
  using TricubicBsplineTraits<T>::Lattice;
  std::vector<PosType> Centers;

  /** default constructure
   *
   * Set Rcut2 to a large number so that everything counts
   */
  TricubicBsplineSPOSet()
  {
    Rcut2=1e6;
  }

  ~TricubicBsplineSPOSet()
  { }

  inline void setTwistAngle(const PosType& tangle)
  { }

  inline void setGrid(const GridType& knots)
  {
    bKnots=knots;
  }

  ///empty reset
  void resetParameters(VarRegistry<real_type>& vlist)
  {
    ///DO NOTHING FOR NOW
  }

  inline void setGrid(real_type xi, real_type xf,
                      real_type yi, real_type yf, real_type zi, real_type zf,
                      int nx, int ny, int nz,
                      bool interp=true, bool periodic=true,bool openend=true)
  {
    bKnots.setGrid(0.0,1.0,0.0,1.0,0.0,1.0,nx,ny,nz,interp,periodic,openend);
  }

  /** add a orbital
   * @param i index of the orbital
   * @param data input data
   * @param curP interpolated data
   */
  void add(int i, const PosType& c, const StorageType& data, StorageType* curP)
  {
    bKnots.Init(data,*curP);
    Centers.push_back(c);
    P.push_back(curP);
  }

  void add(int i, const PosType& c, StorageType* curP)
  {
    //Already exists
    if(i<Centers.size())
      return;
    Centers.push_back(c);
    P.push_back(curP);
  }

  template<typename PV>
  inline void evaluate(const PosType& r, PV& vals)
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.Find(ru[0],ru[1],ru[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
        vals[j]=0.0;//numeric_limits<T>::epsilon();
      else
        vals[j]=bKnots.evaluate(*P[j]);
    }
  }

  template<typename PV, typename GV>
  inline void
  evaluate(const PosType& r, PV& vals, GV& grads, PV& laps)
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    TinyVector<value_type,3> gu;
    Tensor<value_type,3> hess;
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        vals[j]=0.0;//numeric_limits<T>::epsilon();
        grads[j]=0.0;
        laps[j]=0.0;
      }
      else
      {
        vals[j]=bKnots.evaluate(*P[j],gu,hess);
        grads[j]=dot(Lattice.G,gu);
        laps[j]=trace(hess,GGt);
        //vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
      }
    }
  }

  template<typename PM, typename GM>
  inline void
  evaluate(const PosType& r, int i, PM& vals, GM& grads, PM& laps)
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    TinyVector<value_type,3> gu;
    Tensor<value_type,3> hess;
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        vals(j,i)=0.0; //numeric_limits<T>::epsilon();
        grads(i,j)=0.0;
        laps(i,j)=0.0;
      }
      else
      {
        vals(j,i)=bKnots.evaluate(*P[j],gu,hess);
        grads(i,j)=dot(Lattice.G,gu);
        laps(i,j)=trace(hess,GGt);
        //vals(j,i)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
      }
    }
  }

private:
  GridType bKnots;
  std::vector<const StorageType*> P;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2013 $   $Date: 2007-05-22 16:47:09 -0500 (Tue, 22 May 2007) $
 * $Id: TricubicBsplineSet.h 2013 2007-05-22 21:47:09Z jnkim $
 ***************************************************************************/
