#include "QMCTools/QMCFiniteSize/QMCFiniteSize.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include <iostream>
#include <cmath>
#include "Configuration.h"
#include "einspline/bspline_eval_d.h"
#include "einspline/nubspline_eval_d.h"
#include "einspline/nugrid.h"
#include "einspline/nubspline_create.h"
namespace qmcplusplus
{

QMCFiniteSize::QMCFiniteSize():
skparser(NULL),ptclPool(NULL),myRcut(0.0),myConst(0.0),P(NULL),h(0.0)
{
  //std::cout<<" Whoaa constructing QMCFiniteSize\n";
  IndexType mtheta=10;
  IndexType mphi=10;
  app_log()<<"Building spherical grid. n_theta x n_phi = "<<mtheta<<" x "<<mphi<<endl;
  build_spherical_grid(mtheta,mphi);
  h=0.00001; 
}

QMCFiniteSize::QMCFiniteSize(SkParserBase* skparser_i):
skparser(skparser_i),ptclPool(NULL),myRcut(0.0),myConst(0.0),P(NULL),h(0.0),sphericalgrid(0)
{
  //std::cout<<" Whoaa constructing QMCFiniteSize\n";
  mtheta=80;
  mphi=80;
  h=0.00001;
   

}

bool QMCFiniteSize::validateXML()
{
  xmlXPathContextPtr m_context = XmlDocStack.top()->getXPathContext();
  xmlNodePtr cur = XmlDocStack.top()->getRoot()->children;

  while(cur!=NULL)
  {
    std::string cname((const char*)cur->name);
    bool inputnode=true;
    if(cname == "particleset")
    {
      ptclPool.put(cur);
    }
    else if(cname == "wavefunction")
    {
      wfnPut(cur);
    }
    else if(cname == "include")
    {
      //file is provided
      const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"href");
      if(a)
      {
        pushDocument((const char*)a);
        inputnode = processPWH(XmlDocStack.top()->getRoot());
        popDocument();
      }
    }
    else if(cname == "qmcsystem")
    {
      processPWH(cur);
    }
    else
    {
    }
   // if(inputnode)
   //   lastInputNode=cur;
    cur=cur->next;
  }
  
   
  app_log() << "=========================================================\n";
  app_log() << " Summary of QMC systems \n";
  app_log() << "=========================================================\n";
  ptclPool.get(app_log());
}
bool QMCFiniteSize::wfnPut(xmlNodePtr cur)
{
  std::string id("psi0"), target("e"), role("extra");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(id,"id");
  pAttrib.add(id,"name");
  pAttrib.add(target,"target");
  pAttrib.add(target,"ref");
  pAttrib.add(role,"role");
  pAttrib.put(cur);
  ParticleSet* qp = ptclPool.getParticleSet(target);
   
  {//check ESHDF should be used to initialize both target and associated ionic system
    xmlNodePtr tcur=cur->children;
    while(tcur != NULL)
    { //check <determinantset/> or <sposet_builder/> to extract the ionic and electronic structure
      std::string cname((const char*)tcur->name);
      if(cname == OrbitalBuilderBase::detset_tag || cname =="sposet_builder")
      { 
        qp=ptclPool.createESParticleSet(tcur,target,qp);
      }
      tcur=tcur->next;
    }
  }
  
  //std::cout << "=========================================================\n";
 // std::cout << " Summary of Unit Cell Geometry \n";
 // std::cout << "=========================================================\n";
 // qp->Lattice.print(std::cout);

}

void QMCFiniteSize::getSkInfo(UBspline_3d_d* spline, vector<RealType>& symmatelem)
{
  symmatelem.resize(6);
  RealType sx(0), sy(0),sz(0), sxy(0),sxz(0),syz(0);
  RealType h2=h*h;

  PosType disp;
  PosType disp_lat;
  
  disp[0]=h;disp[1]=0;disp[2]=0;
  disp_lat=P->Lattice.k_unit(disp);
  eval_UBspline_3d_d(spline,disp_lat[0],disp_lat[1],disp_lat[2],&sx);
  
  disp[0]=0;disp[1]=h;disp[2]=0;
  disp_lat=P->Lattice.k_unit(disp);
  eval_UBspline_3d_d(spline,disp_lat[0],disp_lat[1],disp_lat[2],&sy);
  
  disp[0]=0;disp[1]=0;disp[2]=h;
  disp_lat=P->Lattice.k_unit(disp);
  eval_UBspline_3d_d(spline,disp_lat[0],disp_lat[1],disp_lat[2],&sz);
  
  disp[0]=h;disp[1]=h;disp[2]=0;
  disp_lat=P->Lattice.k_unit(disp);
  eval_UBspline_3d_d(spline,disp_lat[0],disp_lat[1],disp_lat[2],&sxy);
  
  disp[0]=h;disp[1]=0;disp[2]=h;
  disp_lat=P->Lattice.k_unit(disp);
  eval_UBspline_3d_d(spline,disp_lat[0],disp_lat[1],disp_lat[2],&sxz);
  
  disp[0]=0;disp[1]=h;disp[2]=h;
  disp_lat=P->Lattice.k_unit(disp);
  eval_UBspline_3d_d(spline,disp_lat[0],disp_lat[1],disp_lat[2],&syz);
  
  symmatelem[0]=sx/h2;
  symmatelem[1]=sy/h2;
  symmatelem[2]=sz/h2;
  symmatelem[3]=0.5*(sxy-sx-sy)/h2;
  symmatelem[4]=0.5*(sxz-sx-sz)/h2;
  symmatelem[5]=0.5*(syz-sy-sz)/h2;

}


RealType QMCFiniteSize::sphericalAvgSk(UBspline_3d_d* spline, RealType k)
{
  RealType sum=0.0;
  RealType val=0.0;
  PosType kvec(0);
  IndexType ngrid=sphericalgrid.size();
//  IndexType numpoints=0;
  for(IndexType i=0; i<ngrid; i++)
  {
    
    kvec=k*sphericalgrid[i];
//    if( kvec[0]>gridx.upper_bound || kvec[0]<gridx.lower_bound)
//     continue;
//    if( kvec[1]>gridy.upper_bound || kvec[1]<gridy.lower_bound)
//      continue;
//    if( kvec[2]>gridz.upper_bound || kvec[2]<gridz.lower_bound)
//      continue;
    eval_UBspline_3d_d(spline,kvec[0],kvec[1],kvec[2],&val);
    sum+=val;
//    numpoints++;
  }

  return sum/RealType(ngrid);
}
bool QMCFiniteSize::processPWH(xmlNodePtr cur)
{
  //return true and will be ignored
  if(cur == NULL)
    return true;
  bool inputnode=true;
  //save the root to grep @tilematrix
  xmlNodePtr cur_root=cur;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "simulationcell")
    {
      ptclPool.putLattice(cur);
    }
    else if(cname == "particleset")
    {
      ptclPool.putTileMatrix(cur_root);
      ptclPool.put(cur);
    }
    else if(cname == "wavefunction")
    {
      wfnPut(cur);
    }
//    else if(cname == "hamiltonian")
//    {
//      hamPool->put(cur);
//    }
    else
      //add to m_qmcaction
    {
//      inputnode=false;
//      m_qmcaction.push_back(std::pair<xmlNodePtr,bool>(xmlCopyNode(cur,1),false));
    }
    cur=cur->next;
  }
  //flush
//  app_log().flush();
  return inputnode;
}
void QMCFiniteSize::initBreakup()
{
  
  app_log() << "=========================================================\n";
  app_log() << " Initializing Long Range Breakup (Esler) \n";
  app_log() << "=========================================================\n";
  P =ptclPool.getParticleSet("e");
  AA = LRCoulombSingleton::getHandler(*P);
  //AA->initBreakup(*PtclRef);
  //myConst=evalConsts();
  myRcut=AA->get_rc();//Basis.get_rc();
  if(rVs==0)
  {
    rVs = LRCoulombSingleton::createSpline4RbyVs(AA,myRcut,myGrid);
  }
}

UBspline_3d_d* QMCFiniteSize::getSkSpline(RealType limit)
{
  KContainer Klist = P->SK->KLists;
  
  vector<TinyVector<int,OHMMS_DIM> > kpts = Klist.kpts;

  vector<RealType> sk(kpts.size()), skerr(kpts.size());

  skparser->get_sk(sk,skerr);
  if(skparser->is_normalized()==false)
  {
    IndexType Ne=P->getTotalNum();
    for(IndexType i=0; i<sk.size(); i++)
    {
      sk[i]/=RealType(Ne);
      skerr[i]/=RealType(Ne);
    }
  }
 //Grid_t gridx;
 // Grid_t gridy;
 // Grid_t gridz;

  skparser->get_grid(gridx,gridy,gridz);

  Ugrid esgridx;
  Ugrid esgridy;
  Ugrid esgridz;

  
  //get the einspline grids.
  esgridx = gridx.einspline_grid();
  esgridy = gridy.einspline_grid();
  esgridz = gridz.einspline_grid();
  
  //setup the einspline boundary conditions.  
  BCtype_d bcx;
  BCtype_d bcy;
  BCtype_d bcz;

  
  //This piece iterates through S(k) and sets 
  //pieces beyond the k-cutoff equal to 1.  
  //A violent approximation if S(k) is not converged, but
  //better than S(k)=0.   
  double kc=AA->get_kc();
  double kcutsq=kc*kc;
 
  for(int i=int(gridx.lower_bound), skindex=0; i<=int(gridx.upper_bound); i++)
    for(int j=int(gridy.lower_bound); j<=int(gridy.upper_bound); j++)
      for(int k=int(gridz.lower_bound); k<=int(gridz.upper_bound); k++)
      {
        PosType v;
        v[0]=i;
        v[1]=j;
        v[2]=k;
        RealType ksq=P->Lattice.ksq(v);
                
        if(ksq>kcutsq) sk[skindex]=limit;
        skindex++;
      }
  //No particular BC's on the edge of S(k).  
   
  bcx.lCode=NATURAL;
  bcx.rCode=NATURAL;
  bcx.lVal=1.0;
  bcx.rVal=1.0;

  bcy.lCode=NATURAL;
  bcy.rCode=NATURAL;
  bcy.lVal=1.0;
  bcy.rVal=1.0;
  
  bcz.lCode=NATURAL;
  bcz.rCode=NATURAL;
  bcz.lVal=1.0;
  bcz.rVal=1.0;

  UBspline_3d_d * spline = create_UBspline_3d_d(esgridx,esgridy,esgridz,
                                                bcx, bcy, bcz,
                                                sk.data());

  return spline;
  
}

void QMCFiniteSize::build_spherical_grid(IndexType mtheta, IndexType mphi)
{
  sphericalgrid.resize(mtheta*mphi);
  RealType dphi=2*M_PI/RealType(mphi-1);
  RealType dcostheta=2.0/RealType(mtheta-1);
  
  PosType tmp;


  for(IndexType i=0; i<mtheta; i++)
    for(IndexType j=0; j<mphi; j++)
    {
      IndexType gindex=i*mtheta+j;
      
      RealType costheta=-1.0+dcostheta*i;
      RealType theta=std::acos(costheta);
      RealType sintheta=std::sin(theta);
      
      RealType phi=dphi*j;
      RealType sinphi=std::sin(phi);
      RealType cosphi=std::cos(phi);
      tmp[0]=sintheta*cosphi;
      tmp[1]=sintheta*sinphi;
      tmp[2]=costheta;
      
      //we do this last transformation because S(k) is splined in the lattice 
      //kvector units, and not the cartesian ones.  
      sphericalgrid[gindex]=P->Lattice.k_unit(tmp);

    }

  
}

NUBspline_1d_d* QMCFiniteSize::spline_clamped(vector<RealType>& grid, vector<RealType>& vals, RealType lVal, RealType rVal)
{
  NUgrid* grid1d = create_general_grid(grid.data(),grid.size());

  BCtype_d xBC;
  xBC.lVal=lVal;
  xBC.rVal=rVal;
  xBC.lCode=DERIV1;
  xBC.rCode=DERIV1;
  
  return create_NUBspline_1d_d(grid1d,xBC,vals.data());
}

//Integrate the spline using Simpson's rule.  For Bsplines, this should be exact
//provided your delta is smaller than the smallest bspline mesh spacing.
RealType QMCFiniteSize::integrate_spline(NUBspline_1d_d* spline, RealType a, RealType b, IndexType N)
{
  N=N+3-N%3;
  IndexType ninterv=N/3;

  RealType del=(b-a)/RealType(N);
  RealType s=0.0;
  RealType val=0.0;

  for (IndexType i=0; i<ninterv-1; i++)
  {
    RealType x0(0),x1(0),x2(0),x3(0);
    x0=a+3*i*del;
    x1=a+3*i*del + 1*del;
    x2=a+3*i*del + 2*del;
    x3=a+3*i*del + 3*del;

    
    eval_NUBspline_1d_d(spline,x0,&val);
    s+=val;
 
    eval_NUBspline_1d_d(spline,x1,&val);
    s+=3.0*val;
 
    eval_NUBspline_1d_d(spline,x2,&val);
    s+=3.0*val;
 
    eval_NUBspline_1d_d(spline,x3,&val);
    s+=val;
 
  }
  return s*3.0*del/8.0;
  
  
}

bool QMCFiniteSize::execute()
{
  //Initialize the long range breakup.  For now, do the Esler method.  
  initBreakup();
  KContainer Klist = P->SK->KLists;
  vector<TinyVector<int,OHMMS_DIM> > kpts= Klist.kpts;  //These are in reduced coordinates.
                                     //Easier to spline, but will have to convert 
                                     //for real space integration.  
                                     
  IndexType Ne = P->getTotalNum();
  vector<RealType> sk(kpts.size()),skerr(kpts.size());
 
  if(!skparser->has_grid())
    skparser->set_grid(kpts);

  cout<<"Grid computed\n";
  sk=skparser->get_sk_raw();
  skerr=skparser->get_skerr_raw();
  
  if(skparser->is_normalized()==false)
  {
    for(int i=0; i<sk.size(); i++)
    {
      sk[i]/=RealType(Ne);
      skerr[i]/=RealType(Ne); 
    }
  } 
  //This is the \frac{1}{Omega} \sum_{\mathbf{k}} \frac{v_k}{2} S(\mathbf{k} term.  
  double V=0.5*AA->evaluate_w_sk(Klist.kshell,sk.data());
  vector<RealType> vsk(kpts.size());
  vector<RealType> vsk_1d(Klist.kshell.size());
  
  app_log()<<"#KSHELL_START#\n";
  for(int ks=0; ks<Klist.kshell.size(); ks++)
  {
    app_log() << Klist.kshell[ks] << endl;
  }
  app_log()<<"#KSHELL_STOP#\n";

  for(int ks=0,ki=0; ks<Klist.kshell.size()-1; ks++)
  {
    RealType u=0;
    RealType n=0;
    for(; ki<Klist.kshell[ks+1]; ki++)
    {  
       //u += sk[ki]*AA->Fk[ki];
       u += sk[ki];
       n++;
    }
    
    if(n!=0)
      vsk_1d[ks]=u/n;
    else
      vsk_1d[ks]=0;
  }

  app_log()<<"Spherically Averaged S(k): Raw\n";
  app_log()<<"#SK_RAW_START#\n";
  app_log()<<"k S(k) vk\n";
  for(int ks=0; ks<Klist.kshell.size()-1; ks++)
    app_log()<<std::sqrt(Klist.ksq[Klist.kshell[ks]])<<" "<<vsk_1d[ks]<<" "<<AA->Fk_symm[ks]<<endl;
  app_log()<<"#SK_RAW_STOP#\n";
  UBspline_3d_d* sk3d_spline =  getSkSpline(vsk_1d[Klist.kshell.size()-2]);
//  for (int i=0; i<sk.size(); i++)
//  {
//    PosType mykpt = kpts[i];
//    RealType dval; 
//    eval_UBspline_3d_d(doop,mykpt[0],mykpt[1],mykpt[2],&dval);
//    app_log()<<i<<" "<<sk[i]<<" "<<dval<<endl;
//  }
  
  RealType xstart(-4.0),xend(4.0);
  RealType ystart(-4.0),yend(4.0);
  RealType zstart(-4.0),zend(4.0);
 
  RealType dx,dy,dz;
  RealType nx(16),ny(16),nz(16);

  dx=(xend-xstart)/RealType(nx);
  dy=(yend-ystart)/RealType(ny);
  dz=(zend-zstart)/RealType(nz);

  app_log()<<"Building spherical grid. n_theta x n_phi = "<<mtheta<<" x "<<mphi<<endl;
  build_spherical_grid(mtheta,mphi);
  app_log()<<endl;
  vector<RealType> Amat;
  getSkInfo(sk3d_spline,Amat);
 
  app_log() << "=========================================================\n";
  app_log() << " S(k) Info \n";
  app_log() << "=========================================================\n";
  app_log() << "S(k) anisotropy near k=0\n";
  app_log() << "------------------------\n";
  app_log() << "  a_xx = "<<Amat[0]<<endl;
  app_log() << "  a_yy = "<<Amat[1]<<endl;
  app_log() << "  a_zz = "<<Amat[2]<<endl;
  app_log() << "  a_xy = "<<Amat[3]<<endl;
  app_log() << "  a_xz = "<<Amat[4]<<endl;
  app_log() << "  a_yz = "<<Amat[5]<<endl;
  app_log() << "------------------------\n";


  RealType b = (Amat[0]+Amat[1]+Amat[2])/3.0;
  
  app_log() << "Spherically averaged S(k) near k=0\n";
  app_log() << "S(k)=b*k^2   b = "<<b<<endl;
  app_log() << "------------------------\n";
  app_log()<<endl;
 

  app_log()<<"Now to spherically average with the spline\n";
  RealType kmax=AA->get_kc();
  RealType nk=100;
  RealType kdel=kmax/(nk-1.0);

  // output splined S(k)
  app_log()<< "Spherically Averaged S(k):  Splined\n";
  app_log()<< "#SK_SPLINE_START#\n";
  app_log()<< "k S(k)\n";
  for(int k=0;k<nk; k++)
  {
    RealType kval=kdel*k;
    app_log()<<kval<<" "<<sphericalAvgSk(sk3d_spline,kval)<<endl;
  }
  app_log()<< "#SK_SPLINE_STOP#\n";

  // output vlr(k)
  IndexType ngrid = Klist.kshell.size();
  app_log()<< "#VK_START#\n";
  app_log()<< "k vlr(k)\n";
  for(int ks=0;ks<ngrid-1; ks++)
  {
    RealType kval = std::sqrt(Klist.ksq[Klist.kshell[ks]]);
    app_log()<<kval<<" "<<AA->Fk_symm[ks]<<endl;
  }
  app_log()<< "#VK_STOP#\n";
 
  vector<RealType> nonunigrid1d(ngrid);
  vector<RealType> k2vksk(ngrid); 
  nonunigrid1d[0]=0.0;
  k2vksk[0]=0.0;
  
  RealType myV = 0.0;
  int ki=0;
  for(int ks=0; ks<ngrid-1; ks++)
  {
    
    RealType kval = std::sqrt(Klist.ksq[Klist.kshell[ks]]);
    nonunigrid1d[ks+1]=kval;
    RealType skavg=sphericalAvgSk(sk3d_spline,kval);
    RealType vk = AA->Fk_symm[ks];
    for (;ki<Klist.kshell[ks+1];ki++)
    {
      myV += 0.5*vk*sk[ki];
    }
   
    k2vksk[ks+1]=0.5*vk*skavg*kval*kval;
  }
 
  NUBspline_1d_d* integrand = spline_clamped(nonunigrid1d,k2vksk,0.0,0.0);


  app_log()<< "#INTEGRAND_START#\n";
  app_log()<<"Integrand:\n";
  for(int i=0; i<nk; i++)
  {
    RealType kval = kdel*i;
    RealType dval=0.0;
    eval_NUBspline_1d_d(integrand,kval,&dval); 
    app_log()<<kval<<" "<<dval<<endl;
  }
  app_log()<< "#INTEGRAND_STOP#\n";

  
  //Integrate the spline and compute the thermodynamic limit.
  RealType integratedval = integrate_spline(integrand,0.0,kmax,100);
  RealType intnorm = P->Lattice.Volume/2.0/M_PI/M_PI;  //The volume factor here is because 1/Vol is 
                                                       //included in QMCPACK's v_k.  See CoulombFunctor.
                                                       
  app_log() << "  integratedval = " << integratedval << endl;
  app_log() << "  intnorm = " << intnorm << endl;
  app_log() << "  V = " << V << endl;
  app_log() << "  myV = " << myV << endl;

  app_log() << "=========================================================\n";
  app_log() << " Finite Size Correction:  2-bdy Potential \n";
  app_log() << "=========================================================\n";
  app_log() << "  (All units are in mHa/electron)\n";
  
  RealType vfscorr = intnorm*integratedval - V;
  RealType rho = RealType(Ne)/P->Lattice.Volume;
  RealType rs = std::pow(3.0/(4*M_PI)*P->Lattice.Volume/RealType(Ne),1.0/3.0);
 
  app_log() << "--------------------------------------------------------\n";
  app_log() << "System summary:\n";
  app_log() << "    Ne    = "<<Ne<<endl;
  app_log() << "   Vol    = "<<P->Lattice.Volume<<" a0^3"<<endl;
  app_log() << "  Ne/V    = "<<rho<<" 1/a0^3"<<endl;
  app_log() << " rs/a0    = "<<rs<<endl;
  app_log() << "--------------------------------------------------------\n";
  RealType vlo = 2*M_PI*rho*b/RealType(Ne);
  app_log() << "Finite Size Corrections:\n";
  app_log() << "  V_LO   = "<<vlo*1000<<" mHa/electron\n";
  app_log() << " V_INT   = "<<vfscorr*1000<<" mHa/electron\n";
  app_log() << "\n";
  app_log() << " Using leading order S(k) data, a crude estimate for T_LO:\n";
  RealType tlo = 1.0/RealType(Ne*b*8);
  app_log() << "  T_LO   = "<<tlo*1000<<" mHa/electron\n\n";
  app_log() << " **DISCLAIMER** This is a crude estimate, and can sometimes be\n";
  app_log() << " quite off.  Estimating the T correction from a long range jastrow\n";
  app_log() << " or from mixed and pure structure factors is strongly encouraged.\n";
//  eval_NUBspline_1d_d(sk_1d,kval,&dval);
   
  //vsksort=
  
//  for(int i=0; i<sk.size(); i++)
//    cout<<Klist.kpts_cart[i][0]<<" "<<Klist.kpts_cart[i][1]<<" "<<Klist.kpts_cart[i][2]<<" "<<sk[i]<<" "<<skerr[i]<<endl;
   
   
}

}
