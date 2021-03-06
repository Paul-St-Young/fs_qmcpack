//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Configuration.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
void ReportEngine::echo(xmlNodePtr cur, bool recursive)
{
  if(cur==NULL)
    return;
  app_log()<< "<input node=\""<<(const char*)(cur->name) << "\"";
  xmlAttrPtr att = cur->properties;
  char atext[1024];
  while(att != NULL)
  {
    sprintf(atext,"  %s=\"%s\"",(const char*)(att->name),(const char*)(att->children->content));
    app_log() << atext;
    att = att->next;
  }
  app_log() << "/>\n";
}

bool ReportEngine::DoOutput = false;

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2468 $   $Date: 2008-02-22 09:27:30 -0500 (Fri, 22 Feb 2008) $
 * $Id: Communicate.h 2468 2008-02-22 14:27:30Z jnkim $
 ***************************************************************************/
