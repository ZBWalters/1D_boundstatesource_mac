#include "./classes.h"
#include "./globalbasis.h"

cmplx* __restrict globalbasis::boundstate(int nstate, rl& En){
  //calculate the nth bound state for the given potential
  //nstate=0 corresponds to ground state

  char uplo='U';
  int kd=elementbases[0]->order;
  int ldab=kd+1;
  rl* __restrict Hmat=new rl[nbf*nbf];
  rl* __restrict Omat=new rl[nbf*nbf];
  for(int i=0;i<nbf*nbf;i++){
    Hmat[i]=0.;
  }
  for(int i=0;i<nbf*nbf;i++){
    Omat[i]=0.;
  }

  cout << "setting up bandH,bandO\n";
  for(int eltnum=0;eltnum<nelts;eltnum++){
    derivedbasis* tmpbasis=elementbases[eltnum];
    int tmporder=tmpbasis->order;
    cmplx* __restrict tmph=tmpbasis->Hmat;
    cmplx* __restrict tmpo=tmpbasis->omat;
    cmplx* __restrict tmpnablasq=tmpbasis->nablasqmat;
    cmplx* __restrict tmpV=tmpbasis->Vmat;
    //cout << "tmpo\t"<<eltnum<<"\t"<<tmpbasis->xmax<<"\t"<<tmpbasis->xmin<<"\n";
    //printmat(tmpo,tmporder,tmporder);
    //cout << "\n";
    //
    //cout << "tmpnablasq\t"<<eltnum<<"\t"<<tmpbasis->xmax<<"\t"<<tmpbasis->xmin<<"\n";
    //printmat(tmpnablasq,tmporder,tmporder);
    //cout << "\n";
    //for hbandindx, i
    for(int j=0;j<tmporder;j++){
      int indx2=eltfirstindices[eltnum]+j;
      //cout << "indx2\t"<<indx2<<"\n";
      for(int i=0;i<=j;i++){
	int indx1=eltfirstindices[eltnum]+i;
	//cout << "indx1\t"<<indx1<<"\n";
	//cout << "hbindx\t"<<hbandarrayindx(uplo,indx1,indx2,kd,nbf)<<"\n";
	Hmat[arrayindx(indx1,indx2,nbf,nbf)]+=
	  real(-0.5*tmpnablasq[arrayindx(i,j,tmporder,tmporder)]+tmpV[arrayindx(i,j,tmporder,tmporder)]); //tmph[arrayindx(i,j,tmporder,tmporder)];
	Omat[arrayindx(indx1,indx2,nbf,nbf)]+=
	  real(tmpo[arrayindx(i,j,tmporder,tmporder)]);
      }
    }
  }
  //cout << "set up bandO\n";
  //printmat_scipy(bandO,kd+1,nbf);
  //cout << "\n";
  //cout << "set up bandH\n";
  //printmat_scipy(bandH,kd+1,nbf);


  int itype=1;
  char jobz='V';//'V';//calculate eigenvectors
  char range='I';//find eigenvalues in integer range
  int lda=nbf;
  int ldb=nbf;
  rl vl=0.;//lower range of eigenvalue search (not referenced)
  rl vu=0.;//upper range of eigenvalue search
  int il=1;//nstate+1;
  int iu=nstate+1;
  rl abstol=1.e-6;//eigenvalue tolerance
  int nev=iu-il+1;//number of eigenvalues found
  rl* evals=new rl[nbf];
  rl* __restrict evecs=new rl[nbf*nev];
  int ldz=nbf;
  int lwork=-1;
  rl* tmpwork=new rl[1];
  int* iwork=new int[5*nbf];
  int* ifail=new int[nbf];
  int info=0.;

	dsygvx_(&itype, &jobz, &range, &uplo, &nbf, Hmat, &lda, Omat, &ldb, &vl, &vu,
			&il, &iu, &abstol, &nev, evals, evecs, &ldz, tmpwork, &lwork, iwork,
			ifail,  &info);
  lwork=int(tmpwork[0]);
  delete [] tmpwork;
  rl* work=new rl[lwork];
	dsygvx_(&itype, &jobz, &range, &uplo, &nbf, Hmat, &lda, Omat, &ldb, &vl, &vu,
			&il, &iu, &abstol, &nev, evals, evecs, &ldz, work, &lwork, iwork, ifail,
			&info);
  
  
  //copy results into return values
  cmplx* __restrict retarray=new cmplx[nbf];
  for(int i=0;i<nbf;i++){
    retarray[i]=evecs[nbf*(nev-1)+i];//evecs[i];
  }
  En=evals[nstate];


  //delete temporary matrices
  delete [] Hmat;
  delete [] Omat;
  delete [] evals;
  delete [] evecs;
  delete [] work;
  delete [] iwork;
  delete [] ifail;


  //zero function in any interval where the exterior complex scaling is in use
  for(int eltnum=0;eltnum<nelts;eltnum++){
    derivedbasis* tmpbasis=elementbases[eltnum];
    int tmporder=tmpbasis->order;
    if(tmpbasis->ecsflag){
      for(int i=0;i<tmporder;i++){
	retarray[eltfirstindices[eltnum]+i]=0.;
      }
    }
  }
  
  return retarray;
}

cmplx* __restrict globalbasis::functimespsi(cmplx (*func)(rl x), cmplx* psi_in){
  cmplx* __restrict retvec=initializearray(nbf);
  for(int eltnum=0;eltnum<nelts;eltnum++){
    derivedbasis* tmpbasis=elementbases[eltnum];
    cmplx* tmpfuncpsi=tmpbasis->functimespsi(func,
					     psi_in+eltfirstindices[eltnum]);
    for(int i=0;i<tmpbasis->order;i++){
      retvec[eltfirstindices[eltnum]+i]+=tmpfuncpsi[i];
    }
    delete [] tmpfuncpsi;
  }
  return retvec;
}


////////////////////////////////////////////////////////////////////////
//cmplx* __restrict globalbasis::boundstate(int nstate, rl& En){
//  //calculate the nth bound state for the given potential
//  //nstate=0 corresponds to ground state
//
//  char uplo='U';
//  int kd=elementbases[0]->order;
//  int ldab=kd+1;
//  rl* __restrict bandH=new rl[ldab*nbf];
//  rl* __restrict bandO=new rl[ldab*nbf];
//  for(int i=0;i<ldab*nbf;i++){
//    bandH[i]=0.;
//  }
//  for(int i=0;i<ldab*nbf;i++){
//    bandO[i]=0.;
//  }
//
//  cout << "setting up bandH,bandO\n";
//  for(int eltnum=0;eltnum<nelts;eltnum++){
//    derivedbasis* tmpbasis=elementbases[eltnum];
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
//      int indx2=eltfirstindices[eltnum]+j;
//      //cout << "indx2\t"<<indx2<<"\n";
//      for(int i=0;i<=j;i++){
//	int indx1=eltfirstindices[eltnum]+i;
//	//cout << "indx1\t"<<indx1<<"\n";
//	//cout << "hbindx\t"<<hbandarrayindx(uplo,indx1,indx2,kd,nbf)<<"\n";
//	bandH[hbandarrayindx(uplo,indx1,indx2,kd,nbf)]+=
//	  real(0.5*tmpnablasq[arrayindx(i,j,tmporder,tmporder)]-tmpV[arrayindx(i,j,tmporder,tmporder)]); //tmph[arrayindx(i,j,tmporder,tmporder)];
//	bandO[hbandarrayindx(uplo,indx1,indx2,kd,nbf)]+=
//	  real(tmpo[arrayindx(i,j,tmporder,tmporder)]);
//      }
//    }
//  }
//  //cout << "set up bandO\n";
//  //printmat_scipy(bandO,kd+1,nbf);
//  //cout << "\n";
//  //cout << "set up bandH\n";
//  //printmat_scipy(bandH,kd+1,nbf);
//
//
//  char jobz='N';//'V';
//  char range='I';//find ilth through iuth eigenvalues
//  int ldq=nbf;
//  rl* Q=new rl[ldq*nbf];
//  rl VL=0.;//lower bound of energy range (not used)
//  rl VU=0.;//higher bound of energy range
//  int il=nstate+1;//index of lowest eigenvalue
//  int iu=nstate+1;//index of highest eigenvalue
//  rl abstol=1.e-6;
//  int nev=iu-il+1;//number of eigenvalues
//  rl* evals=new rl[nbf];//eigenvalue array
//  rl* evecs=new rl[nbf*nev];
//  int ldz=nbf;
//  rl* work=new rl[7*nbf];
//  int* iwork=new int[5*nbf];
//  int* ifail=new int[nbf];
//  int info=0.;
//  
//  DSBEVX(&jobz,&range,&uplo,&nbf,&kd,bandH,&ldab,Q,&ldq,&VL,&VU,&il,&iu,
//	 &abstol,&nev,evals,evecs,&ldz,work,iwork,ifail,&info);
//  if(info!=0){
//    cout << "problem with DSBEVX in globalbasis::boundstate!\n";
//  }
//
//  //copy results into return values
//  En=evals[nstate];
//  cmplx* __restrict retarray=new cmplx[nbf];
//  for(int i=0;i<nbf;i++){
//    retarray[i]=evals[i];
//  }
//
//  //delete temporary arrays
//  delete [] bandH;
//  delete [] bandO;
//  delete [] Q;
//  delete [] evals;
//  delete [] evecs;
//  delete [] work;
//  delete [] iwork;
//  delete [] ifail;
//
//
//  
//  return retarray;
//}


//////////////////////////////////////////////////////////////////////////
//cmplx* __restrict globalbasis::boundstate(int nstate, rl& En){
//  //calculate the nth bound state for the given potential
//  //nstate=0 corresponds to ground state
//
//  char uplo='U';
//  int kd=elementbases[0]->order;
//  int ldab=kd+1;
//  rl* __restrict bandH=new rl[ldab*nbf];
//  rl* __restrict bandO=new rl[ldab*nbf];
//  for(int i=0;i<ldab*nbf;i++){
//    bandH[i]=0.;
//  }
//  for(int i=0;i<ldab*nbf;i++){
//    bandO[i]=0.;
//  }
//
//  cout << "setting up bandH,bandO\n";
//  for(int eltnum=0;eltnum<nelts;eltnum++){
//    derivedbasis* tmpbasis=elementbases[eltnum];
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
//      int indx2=eltfirstindices[eltnum]+j;
//      //cout << "indx2\t"<<indx2<<"\n";
//      for(int i=0;i<=j;i++){
//	int indx1=eltfirstindices[eltnum]+i;
//	//cout << "indx1\t"<<indx1<<"\n";
//	//cout << "hbindx\t"<<hbandarrayindx(uplo,indx1,indx2,kd,nbf)<<"\n";
//	bandH[hbandarrayindx(uplo,indx1,indx2,kd,nbf)]+=
//	  real(0.5*tmpnablasq[arrayindx(i,j,tmporder,tmporder)]-tmpV[arrayindx(i,j,tmporder,tmporder)]); //tmph[arrayindx(i,j,tmporder,tmporder)];
//	bandO[hbandarrayindx(uplo,indx1,indx2,kd,nbf)]+=
//	  real(tmpo[arrayindx(i,j,tmporder,tmporder)]);
//      }
//    }
//  }
//  //cout << "set up bandO\n";
//  //printmat_scipy(bandO,kd+1,nbf);
//  //cout << "\n";
//  //cout << "set up bandH\n";
//  //printmat_scipy(bandH,kd+1,nbf);
//
//
//  char jobz='N';//'V';//'V';
//  char range='I';//find ilth through iuth eigenvalues
//  int ldq=nbf;
//  rl* Q=new rl[ldq*nbf];
//  int ldz=nbf;
//  rl* evecs=new rl[ldz*nbf];
//  rl VL=0.;//lower bound of energy range (not used)
//  rl VU=0.;//higher bound of energy range
//  int ka=kd;
//  int kb=kd;
//
//  int il=nstate+1;//index of lowest eigenvalue
//  int iu=nstate+1;//index of highest eigenvalue
//  rl abstol=1.e-6;
//  int nev=iu-il+1;//number of eigenvalues
//  rl* evals=new rl[nbf];//eigenvalue array
//  //rl* evecs=new rl[nbf*nev];
//  rl* work=new rl[7*nbf];
//  int* iwork=new int[5*nbf];
//  int* ifail=new int[nbf];
//  int info=0.;
//  
////  DSBEVX(&jobz,&range,&uplo,&nbf,&kd,bandH,&ldab,Q,&ldq,&VL,&VU,&il,&iu,
////	 &abstol,&nev,evals,evecs,&ldz,work,iwork,ifail,&info);
//  DSBGVX (&jobz,&range,&uplo,&nbf,&ka,&kb,bandH,&ldab,bandO,&ldab,Q,&ldq,
//	  &VL,&VU,&il,&iu,&abstol,&nev,evals,evecs,&ldz,work,iwork,ifail,&info);
//
//  if(info!=0){
//    cout << "problem with DSBGVX in globalbasis::boundstate!\n";
//  }
//
//  //copy results into return values
//  En=evals[nstate];
//  cmplx* __restrict retarray=new cmplx[nbf];
//  for(int i=0;i<nbf;i++){
//    retarray[i]=evals[arrayindx(i,nstate,nbf,1)];
//  }
//
//  //delete temporary arrays
//  delete [] bandH;
//  delete [] bandO;
//  delete [] Q;
//  delete [] evals;
//  delete [] evecs;
//  delete [] work;
//  delete [] iwork;
//  delete [] ifail;
//
//
//  
//  return retarray;
//}
