//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file NewTimer.cpp
 * @brief Implements TimerManager
 */
#include "Utilities/NewTimer.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include <map>
#include <limits>
#include <cstdio>
namespace qmcplusplus
{
TimerManagerClass TimerManager;

void TimerManagerClass::reset()
{
  for (int i=0; i<TimerList.size(); i++)
    TimerList[i]->reset();
}

void
TimerManagerClass::print(Communicate* comm)
{
#if !defined(DISABLE_TIMER)
  std::map<std::string,int> nameList;
  std::vector<double> timeList;
  std::vector<long>   callList;
  //std::vector<int>    callers;
  //timeList.reserve(TimerList.size());
  //callList.reserve(TimerList.size());
  for(int i=0; i<TimerList.size(); ++i)
  {
    NewTimer &timer = *TimerList[i];
    std::map<std::string,int>::iterator it(nameList.find(timer.get_name()));
    if(it == nameList.end())
    {
      int ind=nameList.size();
      nameList[timer.get_name()]=ind;
      timeList.push_back(timer.get_total());
      callList.push_back(timer.get_num_calls());
    }
    else
    {
      int ind=(*it).second;
      timeList[ind]+=timer.get_total();
      callList[ind]+=timer.get_num_calls();
    }
  }
  comm->allreduce(timeList);
  comm->allreduce(callList);
  if(comm->rank() == 0)
  {
    #pragma omp master
    {
      std::map<std::string,int>::iterator it(nameList.begin()), it_end(nameList.end());
      while(it != it_end)
      {
        int i=(*it).second;
        //if(callList[i]) //skip zeros
        fprintf (stderr, "%-40s  %9.4f  %13ld  %16.9f  %12.6f TIMER\n"
        , (*it).first.c_str()
        , timeList[i], callList[i]
        , timeList[i]/(static_cast<double>(callList[i])+numeric_limits<double>::epsilon())
        , timeList[i]/static_cast<double>(omp_get_max_threads()*comm->size()));
        ++it;
      }
    }
  }
#endif
}
}
