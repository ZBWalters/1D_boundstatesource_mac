#include "./classes.h"
#include "./wf.h"

rl wf::psinorm(){
  return gbas->norm(psi);
}

cmplx* __restrict wf::boundstate(int n, rl& En){
  return gbas->boundstate(n,En);
}

void wf::boundstatesetup(int nstate, rl&En, rl t){
  cmplx II=cmplx(0.,1.);
  psi=boundstate(nstate,En);
  //wffile1.open(("boundstate.dat").c_str());
  ofstream wffile1("boundstate.dat");
  wffile1 << this->printstr();
  wffile1.close();
  wffile1.clear();
  arrayscale(psi,nbf,exp(-II*En*t));
}

void wf::boundstatesetup(int nstate, rl& En){
  if(psi!=0){
    delete [] psi;
  }
  psi=boundstate(nstate,En);
  //wffile1.open(("boundstate.dat").c_str());
  ofstream wffile1("boundstate.dat");
  wffile1 << this->printstr();
  wffile1.close();
  wffile1.clear();
}

void wf::boundstatesetup_source(int nstate, rl& En){
  //since the bound state will be operating as a source term, it will
  //be convenient to calculate x*psi rather than psi
  if(psi!=0){
    delete [] psi;
  }
  cmplx* psitmp=boundstate(nstate,En);
  //wffile1.open(("boundstate.dat").c_str());
  ofstream wffile1("boundstate.dat");
  wffile1 << this->printstr();
  wffile1.close();
  wffile1.clear();

  psi=gbas->functimespsi(X_func,psitmp);
  delete [] psitmp;
}



//cmplx* __restrict wf::boundstate(int n, rl& En){
//  
//
////  int kl=gbas->xorder*2;
////  int ku=gbas->xorder*2;
////  int width=2*kl+ku+1;
////  cmplx* __restrict fullbandO=initializearray(nbf*width);
////  cmplx* __restrict fullbandH=initializearray(nbf*width);
////  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
////    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
////    int tmporder=tmpbasis->order;
////    cmplx* __restrict tmpo=tmpbasis->omat;
////    //instead of filling in matrix item by item, do it with blas calls
////    for(int j=0;j<tmporder;j++){
////      cmplx coeff=1.;
////      int indx1=gbas->eltfirstindices[eltnum];
////      int indx2=gbas->eltfirstindices[eltnum]+j;
////      int inc=1;
////      cblas_zaxpy(tmporder,&coeff,tmpo+arrayindx(0,j,tmporder,tmporder),
////                  inc,fullbandO+bandarrayindx(indx1,indx2,ku,kl,nbf),inc);
////    }
////  }
////  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
////    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
////    int tmporder=tmpbasis->order;
////    cmplx* __restrict tmph=tmpbasis->Hmat;
////    //instead of filling in matrix item by item, do it with blas calls
////    for(int j=0;j<tmporder;j++){
////      cmplx coeff=1.;
////      int indx1=gbas->eltfirstindices[eltnum];
////      int indx2=gbas->eltfirstindices[eltnum]+j;
////      int inc=1;
////      cblas_zaxpy(tmporder,&coeff,tmph+arrayindx(0,j,tmporder,tmporder),
////                  inc,fullbandH+bandarrayindx(indx1,indx2,ku,kl,nbf),inc);
////    }
////  }
////
////  //copy fullbandO, fullbandH into band0, bandH
////  int ka=ku;
////  int ldab=ka+1;
////  char uplo='u';
////  cmplx* __restrict bandH=initializearray(ldab*nbf);
////  cmplx* __restrict bandO=initializearray(ldab*nbf);
////  for(int j=0;j<nbf;j++){
////    for(int i=max(0,j-ka);i<=j;i++){
////      bandO[hbandarrayindx(uplo,i,j,ku,nbf)]=
////	fullbandO[bandarrayindx(i,j,ku,kl,nbf)];
////    }
////  }
////  for(int j=0;j<nbf;j++){
////    for(int i=max(0,j-ka);i<=j;i++){
////      bandH[hbandarrayindx(uplo,i,j,ku,nbf)]=
////	fullbandH[bandarrayindx(i,j,ku,kl,nbf)];
////    }
////  }
//
//  
//
//  //copy hmat into bandH, olap into bandO
//
//  char uplo='U';
//  int ka=gbas->elementbases[0]->order;//number of superdiagonals
//  int ldab=ka+1;
//  cout << "ka,nbf\t"<<ka<<"\t"<<nbf<<"\n";
//  cout << "nbf\t"<<nbf<<"\n";
//  cout << "ldab*nbf\t"<<ldab*nbf<<"\n";
//  cmplx* __restrict bandH=initializearray(ldab*nbf);
//  cmplx* __restrict bandO=initializearray(ldab*nbf);
//
//
//
//  cout << "setting up bandH,bandO\n";
//  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
//    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
//    int tmporder=tmpbasis->order;
//    cmplx* __restrict tmph=tmpbasis->Hmat;
//    cmplx* __restrict tmpo=tmpbasis->omat;
//    cmplx* __restrict tmpnablasq=tmpbasis->nablasqmat;
//    cmplx* __restrict tmpV=tmpbasis->Vmat;
//    cout << "tmpo\t"<<eltnum<<"\t"<<tmpbasis->xmax<<"\t"<<tmpbasis->xmin<<"\n";
//    printmat(tmpo,tmporder,tmporder);
//    cout << "\n";
//
//    cout << "tmpnablasq\t"<<eltnum<<"\t"<<tmpbasis->xmax<<"\t"<<tmpbasis->xmin<<"\n";
//    printmat(tmpnablasq,tmporder,tmporder);
//    cout << "\n";
//    //for hbandindx, i
//    for(int j=0;j<tmporder;j++){
//      int indx2=gbas->eltfirstindices[eltnum]+j;
//      //cout << "indx2\t"<<indx2<<"\n";
//      for(int i=0;i<=j;i++){
//	int indx1=gbas->eltfirstindices[eltnum]+i;
//	//cout << "indx1\t"<<indx1<<"\n";
//	//cout << "hbindx\t"<<hbandarrayindx(uplo,indx1,indx2,ka,nbf)<<"\n";
//	bandH[hbandarrayindx(uplo,indx1,indx2,ka,nbf)]+=
//	  real(0.5*tmpnablasq[arrayindx(i,j,tmporder,tmporder)]-tmpV[arrayindx(i,j,tmporder,tmporder)]); //tmph[arrayindx(i,j,tmporder,tmporder)];
//	bandO[hbandarrayindx(uplo,indx1,indx2,ka,nbf)]+=
//	  real(tmpo[arrayindx(i,j,tmporder,tmporder)]);
//      }
//    }
//  }
//  cout << "set up bandO\n";
//  printmat_scipy(bandO,ka+1,nbf);
//  cout << "\n";
//  cout << "set up bandH\n";
//  printmat_scipy(bandH,ka+1,nbf);
//
//  //set up call to zhbgvx for eigenvalues
//  char jobz='N';//'V';//compute eigenvectors
//  char range='I';//compute ILth through IUth evals, evecs
//  int kb=ka;//number of superdiagonals
//  int ldbb=ldab;
//  int ldq=nbf;//number of eigenvectors
//  cmplx* Q=initializearray(ldq*nbf);
//  cout << "initialized Q\n";
//  rl vl=0.01;//lower end of eigenvalue interval (not referenced)
//  rl vu=0.01;//upper end of eigenvalue interval (not referenced)
//  int il=1;//index of smallest eigenvalue to be returned
//  int iu=1;//index of largest eigenvalue to be returned
//  rl abstol=1.e-3;//absolute tolerance for eigenvalues
//  int m=iu-il+1;//total number of eigenvalues found
//  rl*  w=new rl[nbf];//eigenvalues in ascending order on return
//  for(int i=0;i<nbf;i++){
//    w[i]=0.;
//  }
//  cout << "initialized w\n";
//  int ldz=nbf;//leading dimension of array Z (number of eigenvectors)
//  cout << "ldz\t"<<ldz<<"\n";
//  cmplx* Zmat=new cmplx[ldz*nbf];
//  cout <<"initialized Z\n";
//  cmplx* work=new cmplx[nbf];
//  cout << "initialized Z,work\n";
//  rl* rwork=new rl[7*nbf];
//  int* iwork=new int[5*nbf];
//  int* ifail=new int[nbf];
//  int info=0;
//
//  cout << "inputs\n";
//  cout << jobz<<"\t"<<range<<"\t"<<uplo<<"\t"<<nbf<<"\t"<<ka<<"\t"<<kb<<"\n";
//  cout << ldab<<"\t"<<ldbb<<"\t"<<ldq<<"\t"<<vl<<"\t"<<vu<<"\t"<<il<<"\t"<<iu<<"\t"<<abstol<<"\t"<<m<<"\t"<<ldz<<"\n";
////  zhbgvx(&jobz,&range,&uplo,&nbf,&ka,&kb,bandH,&ldab,bandO,&ldbb,Q,&ldq,&vl,&vu,
////	 &il,&iu,&abstol,&m,w,Zmat,&ldz,work,rwork,iwork,ifail,&info);
//  ZHBGVX(&jobz,&range,&uplo,&nbf,&ka,&kb,bandH,&ldab,bandO,&ldbb,Q,&ldq,
//	 &vl,&vu,&il,&iu,&abstol,&m,w,Zmat,&ldz,work,rwork,
//	 iwork,ifail,&info);
//  if(info!=0){
//    cout << "problem with zhbgvx in boundstate!\t"<<info<<"\n";
//  }
//  cout << "back from zhgbvx\n";
//
//  cout <<"eigenenergies\n";
//  
//			
//
//  En=w[n];
//  cout << "En\t"<<En<<"\n";
//  
//  delete [] w;
//  delete [] work;
//  delete [] iwork;
//  delete [] ifail;
//
//  return Zmat;//eigenvector
//
//}
