#include "./classes.h"
#include "./wf.h"

void wfxt::minimizeaction(cmplx* psi0){
  delete [] psixt;
  wfxt_psi0xF0(psi0);
  minimizeaction();
}

void wfxt::wfxt_psi0xF0(cmplx* psi0){
  cout << "wfxt_psi0xF0 "<<nbf<<", "<<torder<<"\n";
  if(psixt!=0){
    delete [] psixt;
  }
  psixt=psi0xF0(nbf,torder,torder,psi0);
}

void wfxt::minimizeaction(){
  //find correction to wfxt which minimizes total action, subject to
  //constraint that psixt(t0) is unaffected
  cmplx* __restrict dpsi0=initializearray(nbf);
  cmplx* __restrict dpsixt=gbas->minimizeactioncorrection(psixt,dpsi0);
  
  cout << "psixt before correction\n";
  printmat(psixt,nbf,torder);

  for(int i=0;i<nbf*torder;i++){
    //cout << "i, psixt[i], dpsixt[i] "<<i<<", "<<psixt[i]<<", "<<dpsixt[i]<<"\n";
    psixt[i]+=dpsixt[i];
    //psixt[i]=dpsixt[i];
  }
  cout << "psixt after correction\n";
  printmat(psixt,nbf,torder);
  delete [] dpsi0;
  delete [] dpsixt;
}

rl wfxt::timestep(rl errortarget){
  rl mindt=1.e-6;
  cmplx* __restrict psilast=initializearray(nbf);
  for(int i=0;i<nbf;i++){
    psilast[i]=psixt[arrayindx(i,torder-1,nbf,torder)];
  }

  cmplx* __restrict psif=psi_f_array();

  cout << "psilast\n";
  printmat(psilast,nbf,1);

  //rl norm=0.;
  //for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
  //  rl eltnorm=gbas->elementbases[eltnum]->norm(&psilast[gbas->eltfirstindices[eltnum]]);
  //  norm+=eltnorm;
  //  //cout << "eltnorm "<<eltnum<<", "<<eltnorm<<"\n";
  //}
  //cout << "norm of last wf component "<<norm<<", "<<errortarget<<"\n";
  //cout << "norm test "<<gbas->norm(psilast)<<"\n";
  rl norm=gbas->norm(psilast);
  rl dnormdn=gbas->innerproduct(psilast,psif);

  rl tmin=gbas->tbasis->tmin;
  rl dt=(gbas->tbasis->tmax - gbas->tbasis->tmin);
  cout << "dt, norm, dnormdn "<<dt<<", "<<norm<<", "<<dnormdn<<"\n";


  rl newdt=dt*pow((pow(errortarget,2)/norm),(1./torder));
  cout << "t, dt, newdt  (norm way)"<<tmin<<", "<<dt<<", "<<newdt<<"\n";

//  newdt=dt*pow(abs(errortarget/dnormdn),(1./(torder-1.)));
//  cout << "t, dt, newdt  (dnormdn way)"<<tmin<<", "<<dt<<", "<<newdt<<"\n";
  

  ////cmplx action=abs(gbas->action(psixt));
  ////cout << "action "<<action<<"\n";
  //newdt=dt*pow((errortarget/abs(action)),(1./(torder-1.)));
  //cout << "dt, newdt  (action way)"<<dt<<", "<<newdt<<"\n";

  delete [] psilast;

  return max(mindt,newdt);  
}
