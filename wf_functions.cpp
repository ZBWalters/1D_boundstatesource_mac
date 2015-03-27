#include "./classes.h"
#include "./wf.h"

void wf::psi0setup(cmplx (*func)(rl x)){
  delete [] psi;
  psi=gbas->psi0setup(func);
}

void wf::psi0setup(cmplx (*func)(rl x, rl kappa, rl rc),rl kappa, rl rc){
  delete [] psi;
  psi=gbas->psi0setup(func,kappa,rc);
}

void wf::print(){
  cout << printstr();
}

void wf::print(int neltpts){
  cout << printstr(neltpts);
}

str wf::printstr(int neltpts){
  //string corresponding to printout of each element evaluated at its
  //integrtion points
  int nelts=gbas->nelts;
  str retstr="";
  for(int i=0;i<nelts;i++){
    retstr+=gbas->elementbases[i]->printstr(neltpts,&psi[eltfirstindices[i]]);
  }
  return retstr;
}

str wf::printstr(){
  //string corresponding to printout of each element evaluated at its
  //integrtion points
  int nelts=gbas->nelts;
  str retstr="";
  for(int i=0;i<nelts;i++){
    retstr+=gbas->elementbases[i]->printstr(&psi[eltfirstindices[i]]);
  }
  return retstr;
}




//////wfxt functions
void wfxt::psi0setup(wf* psi0){
  wfxt_psi0xF0(psi0->psi);
  //delete [] psixt;
  //psixt=initializearray(nbf*torder);
  ////choose psixt to be low order function equal to psi0 at t=0,
  ////0 at t=dt
  //for(int i=0;i<nbf;i++){
  //  psixt[arrayindx(i,0,nbf,torder)]=0.5*psi0->psi[i];
  //  psixt[arrayindx(i,1,nbf,torder)]=-0.5*psi0->psi[i];
  //}
}

wf* wfxt::psi_i(){
  wf* retwf=new wf(gbas,psi_i_array());
  return retwf;
}

cmplx* wfxt::psi_i_array(){
  //return psi at initial time
  rl tmin=gbas->tbasis->tmin;
  cmplx* __restrict tvals=gbas->tbasis->evalbasisfuncs(tmin);

  cmplx* __restrict retwf=initializearray(nbf);
  for(int i=0;i<nbf;i++){
    for(int n=0;n<torder;n++){
      retwf[i]+=psixt[arrayindx(i,n,nbf,torder)]*tvals[n];
    }
  }

  delete [] tvals;
  return retwf;
}

wf* wfxt::psi_f(){
  wf* retwf=new wf(gbas,psi_f_array());
  return retwf;
}

cmplx* wfxt::psi_f_array(){
  //return psi at initial time
  rl tmax=gbas->tbasis->tmax;
  cout << "tmax "<<tmax<<"\n";
  cmplx* __restrict tvals=gbas->tbasis->evalbasisfuncs(tmax);

  //cout << "psixt in psi_f_array\n";
  //printmat(psixt,nbf,torder);

  cout <<"tmin, tmax "<<gbas->tbasis->tmin<<"\t"<<gbas->tbasis->tmax<<"\n";
  cout <<"xmin, xmax "<<gbas->tbasis->xmin<<"\t"<<gbas->tbasis->xmax<<"\n";

  cout << "psif tvals\n";
  printmat(tvals,torder,1);


  //indexing issue: how to reconstruct final wf?
  cmplx* __restrict retwf=initializearray(nbf);
  for(int i=0;i<nbf;i++){
    for(int n=0;n<torder;n++){
      //retwf[i]+=psixt[lsindx(i,n,torder,nbf)]*tvals[n];
      //cout << "i,n,psixt[i,n] "<<i<<", "<<n<<", "<<psixt[arrayindx(i,n,nbf,torder)]<<"\n";
      retwf[i]+=psixt[arrayindx(i,n,nbf,torder)]*tvals[n];
    }
  }
  cout << "retwf\n";
  printmat(retwf,nbf,1);

  delete [] tvals;
  return retwf;
}
