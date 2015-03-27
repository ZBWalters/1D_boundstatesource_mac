#include "./classes.h"
#include "./wf.h"

//The routines in this file project the residual onto generalized eigenvectors
//of Q_R v= lambda U_R v, where Q_R and U_R are the derivative and overlap
//operators in the temporal basis with the zeroth order basis function removed
//(so that Q_R is invertible).  For each temporal eigenvector v, the residuals
//of order higher than 0 are projected onto the eigenvector, and the
//corrections needed to cancel these residuals are calculated and applied.
//The ultimate result is a quickly converging, numerically stable procedure to
//cancel residuals at levels higher than 0.

cmplx* __restrict wf::deltawf_eigenvector(cmplx* __restrict resid, 
					  cmplx Qval,cmplx Uval,
					  cmplx Aval,cmplx Asqval){
  //solve (I Olap alphaQ Q - H alphaU -I p alphaA-Olap AlphaAsq ) dpsi = resid

  cmplx II=cmplx(0.,1.);
  rl ecstheta=gbas->ecstheta;
  int neq=nbf;

  //cout << "inside residualcorrection\n";

  //set up rhsmat to be -resid_{i,n}
  cmplx* __restrict rhsvec=initializearray(neq);
  int inc=1;
  cmplx coeff=1.;
  cblas_zaxpy(neq,&coeff,resid,inc,rhsvec,inc);
  //for(int i=0;i<neq;i++){
  //  rhsvec[i]=-resid[i];
  //}
  //cout << "rhsvec\n";
  //printmat(rhsvec,neq,1);
  //cout << "residxt[n]\n";
  //printmat(residxt[n],neq,1);


  int kl=gbas->xorder;
  int ku=gbas->xorder;
  int width=2*kl+ku+1;

  //copy element hmat, olap into lhsmat
  cmplx* __restrict lhsmat=initializearray(nbf*width);
  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
    int tmporder=tmpbasis->order;
    cmplx* __restrict tmph=tmpbasis->Hmat;
    cmplx* __restrict tmpnablasq=tmpbasis->nablasqmat;
    cmplx* __restrict tmpV=tmpbasis->Vmat;
    cmplx* __restrict tmpo=tmpbasis->omat;
    cmplx* __restrict tmpnabla=tmpbasis->nablamat;
    cmplx* __restrict tmptransnabla=transpose(tmpnabla,tmporder,tmporder);


    //cout << "tmpnabla\n";
    //printmat(tmpnabla,tmporder,tmporder);
    //cout << "tmpnablasq\n";
    //printmat(tmpnablasq,tmporder,tmporder);

//    //set up p.olap, olap.p for use in velocity gauge hamiltonian
    cmplx* __restrict tmponabla=initializearray(tmporder*tmporder);
    cmplx* __restrict tmpnablao=initializearray(tmporder*tmporder);
    cmplx* __restrict tmposq=initializearray(tmporder*tmporder);

    cmplx alpha=1.;
    cmplx beta=0.;
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,tmporder,tmporder,
		tmporder,&alpha,tmpo,tmporder,tmpnabla,tmporder,&beta,
		tmponabla,tmporder);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,tmporder,tmporder,
		tmporder,&alpha,tmpnabla,tmporder,tmpo,tmporder,&beta,
		tmpnablao,tmporder);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,tmporder,tmporder,
		tmporder,&alpha,tmpo,tmporder,tmpo,tmporder,&beta,
		tmposq,tmporder);

    
    //cout << "tmph "<<eltnum<<"\n";
    //printmat(tmph,tmporder,tmporder);
    //cout << "tmpo "<<eltnum<<"\n";
    //printmat(tmpo,tmporder,tmporder);

    //instead of filling in lhsmat item by item, do it with blas calls
    for(int j=0;j<tmporder;j++){
      cmplx* __restrict tmpsum=initializearray(tmporder);
      //cmplx* __restrict tmpsum2=initializearray(tmporder);

      int inc=1;
      cmplx Qcoeff=II*Qval;
      //for arrayindx, i is the fast index, so looping through i can be
      //replaced with blas calls
      cblas_zaxpy(tmporder,&Qcoeff,tmpo+arrayindx(0,j,tmporder,tmporder),
		  inc,tmpsum,inc);
      cmplx Ucoeff2=-Uval;
      cblas_zaxpy(tmporder,&Ucoeff2,tmpV+arrayindx(0,j,tmporder,tmporder),
		  inc,tmpsum,inc);

      cmplx Ucoeff1=.5*Uval;
      cblas_zaxpy(tmporder,&Ucoeff1,tmpnablasq+arrayindx(0,j,tmporder,tmporder),
		  inc,tmpsum,inc);
      cmplx Acoeff=0.5*II*Aval;
      //if(tmpbasis->ecsflag){
      //	Acoeff*=exp(II*ecstheta);
      //}
      cblas_zaxpy(tmporder,&Acoeff,tmpnabla+arrayindx(0,j,tmporder,tmporder),
		  inc,tmpsum,inc);
      Acoeff=-0.5*II*Aval;
      //if(tmpbasis->ecsflag){
      //	Acoeff*=exp(-II*ecstheta);//exp(II*ecstheta);
      //}
      cblas_zaxpy(tmporder,&Acoeff,tmptransnabla+arrayindx(0,j,tmporder,tmporder),
		  inc,tmpsum,inc);
      
      cmplx Asqcoeff=.5*Asqval;
      //if(tmpbasis->ecsflag){
      //	Asqcoeff*=exp(2.*II*ecstheta);
      //}
      cblas_zaxpy(tmporder,&Asqcoeff,tmpo+arrayindx(0,j,tmporder,tmporder),
		  inc,tmpsum,inc);

//      for(int i=0;i<tmporder;i++){	
//	tmpsum2[i]=II*Qval*tmpo[arrayindx(i,j,tmporder,tmporder)]-
//          Uval*(-.5*tmpnablasq[arrayindx(i,j,tmporder,tmporder)]+
//		tmpV[arrayindx(i,j,tmporder,tmporder)])-
//	  Aval*(II*tmpnabla[arrayindx(i,j,tmporder,tmporder)])+ 
//	  .5*Asqval*tmpo[arrayindx(i,j,tmporder,tmporder)];
//      }
//      cout << "tmpsum\n";
//      printmat(tmpsum,tmporder,1);
//      cout << "tmpsum2\n";
//      printmat(tmpsum2,tmporder,1);
//      delete [] tmpsum2;
	
	
      //for bandarrayindx, i is the fast index, so looping through i can be
      //replaced with blas calls
      cmplx sumcoeff=1.;
      int indx1=gbas->eltfirstindices[eltnum];
      int indx2=gbas->eltfirstindices[eltnum]+j;
      int width=2*kl+ku+1;//bandarrayindx decreases with increasing j
      cblas_zaxpy(tmporder,&sumcoeff,tmpsum,inc,
		  lhsmat+bandarrayindx(indx1,indx2,ku,kl,nbf),1);
      delete [] tmpsum;
    }

    //for(int j=0;j<tmporder;j++){
    //  int indx2=gbas->eltfirstindices[eltnum]+j;
    //  for(int i=0;i<tmporder;i++){
    //	int indx1=gbas->eltfirstindices[eltnum]+i; 
    //	cout << "i "<<i<<"\t"<<j<<"\t"<<arrayindx(i,j,tmporder,tmporder)<<"\t"
    //	     <<bandarrayindx(indx1,indx2,ku,kl,nbf)<<"\n";
    //	lhsmat[bandarrayindx(indx1,indx2,ku,kl,nbf)]+=
    //      II*Qval*tmpo[arrayindx(i,j,tmporder,tmporder)]-
    //      Uval*(-.5*tmpnablasq[arrayindx(i,j,tmporder,tmporder)]+
    //		tmpV[arrayindx(i,j,tmporder,tmporder)])-
    //	  Aval*(II*tmpnabla[arrayindx(i,j,tmporder,tmporder)])+
    //	  .5*Asqval*tmpo[arrayindx(i,j,tmporder,tmporder)];
    //  }
    //}

    delete [] tmpnablao;
    delete [] tmponabla;
    delete [] tmposq;
    delete [] tmptransnabla;
  }

  //solve for deltapsin using zgbsv_
  int info=0;
  int* ipiv=new int[neq];
  int nrhs=1;
  int ldab=2*kl+ku+1;

	zgbsv_(&neq, &kl, &ku, &nrhs, reinterpret_cast <Lcmplx*> (lhsmat), &ldab,
			ipiv, reinterpret_cast <Lcmplx*> (rhsvec), &neq, &info);
  if(info!=0){
    cout << "problem with zgbsv_ in residualcorrection! "<<info<<"\n";
  }
  //delete [] AUmat;
  //delete [] AsqUmat;
  delete [] lhsmat;
  delete [] ipiv;

  //test solution
  //cout << "residxt[n]\n";
  //printmat(residxt[n],neq,1);
  //cmplx* Hdeltapsi=gbas->Hpsi(rhsvec);
  //cmplx* Olapdeltapsi=gbas->Olappsi(rhsvec);
  //cmplx* deltaresid=initializearray(neq);
  //for(int i=0;i<neq;i++){
  //  deltaresid[i]=II*Olapdeltapsi[i]*Qval-Hdeltapsi[i]*Uval;
  //}
  //cout << "delta resid\n";
  //printmat(deltaresid,neq,1);
  //end test solution

  return rhsvec;
}

///////////////////////////////////////////////////////////////////////////

//apply corrections in the same way that the linear system is set up
void wf::applydeltapsi(int ncorr,int nmax, cmplx* __restrict deltapsi,
		       cmplx* __restrict * __restrict guessxt,
		       cmplx* __restrict * __restrict residxt, 
		       temporalbasis* tbas ){
  cmplx II=cmplx(0.,1.);
  rl ecstheta=gbas->ecstheta;
  int ot=tbas->order;
  //apply deltapsi vector to entry n of guessxt, adjust residuals by
  //L.deltapsi

  //cout << "inside applydeltapsi\n";
  //cout <<"umat\n";
  //printmat(tbas->umat,ot,ot);
  //cout <<"qmat\n";
  //printmat(tbas->qmat,ot,ot);

  //apply deltapsi vector to entry n of guessxt
  //for(int i=0;i<nbf;i++){
  //  guessxt[ncorr][i]+=deltapsi[i];
  //}
  int inc=1;
  cmplx one=1.;
  cblas_zaxpy(nbf,&one,deltapsi,inc,guessxt[ncorr],inc);


  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
    int tmporder=tmpbasis->order;
    int eltindx=gbas->eltfirstindices[eltnum];
    cmplx* __restrict tmph=tmpbasis->Hmat;
    cmplx* __restrict tmpnablasq=tmpbasis->nablasqmat;
    cmplx* __restrict tmpV=tmpbasis->Vmat;
    cmplx* __restrict tmpo=tmpbasis->omat;
    cmplx* __restrict tmpnabla=tmpbasis->nablamat;
    cmplx* __restrict tmptransnabla=transpose(tmpnabla,tmporder,tmporder);


    //cout << "tmpnabla\n";
    //printmat(tmpnabla,tmporder,tmporder);
    //cout << "tmpnablasq\n";
    //printmat(tmpnablasq,tmporder,tmporder);


    cmplx alpha=1.;
    cmplx beta=0.;

    //find product of Hamiltonian & overlap arrays within temporary basis
    cmplx* __restrict tmpodeltapsi=tmpbasis->Xmatpsi(tmpbasis->omat,
						&deltapsi[eltindx]);
    cmplx* __restrict tmpVdeltapsi=tmpbasis->Xmatpsi(tmpbasis->Vmat,
						&deltapsi[eltindx]);
    cmplx* __restrict tmpnablasqdeltapsi=tmpbasis->Xmatpsi(tmpbasis->nablasqmat,
						&deltapsi[eltindx]);
    cmplx* __restrict tmpnabladeltapsi=tmpbasis->Xmatpsi(tmpbasis->nablamat,
						&deltapsi[eltindx]);

    cmplx* __restrict tmptransposenabladeltapsi=tmpbasis->
      Xmatpsi(tmptransnabla,&deltapsi[eltindx]);
    
    //sum L.psi for current element & apply to residxt
    for(int n=0;n<nmax;n++){
      cmplx* __restrict tmpdeltaresid=initializearray(tmporder);
      cmplx coeff=0.;
      //      if(n<=max(ncorr,2)){
      cmplx qval=tbas->qmat[arrayindx(n,ncorr,ot,ot)];
      coeff=II*qval;
      cblas_zaxpy(tmporder,&coeff,tmpodeltapsi,inc,tmpdeltaresid,inc);
      //}
      //if(n>=(max(0,ncorr-2)) and (n<min(nmax,ncorr+2))){
      cmplx uval=tbas->umat[arrayindx(n,ncorr,ot,ot)];
      coeff=-uval;
      cblas_zaxpy(tmporder,&coeff,tmpVdeltapsi,inc,tmpdeltaresid,inc);

      coeff=0.5*uval;
      cblas_zaxpy(tmporder,&coeff,tmpnablasqdeltapsi,inc,tmpdeltaresid,inc);
	//}

      cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)];
      if(tmpbasis->ecsflag){
	coeff=0.5*II*Aval;//*exp(II*ecstheta);
      }
      else{
	coeff=0.5*II*Aval;
      }
      //coeff*=-1.;
      cblas_zaxpy(tmporder,&coeff,tmpnabladeltapsi,inc,tmpdeltaresid,inc);

      if(tmpbasis->ecsflag){
	coeff=-0.5*II*Aval;//*exp(-II*ecstheta);//exp(II*ecstheta);
      }
      else{
	coeff=-0.5*II*Aval;
      }
      //coeff*=-1.;
      cblas_zaxpy(tmporder,&coeff,tmptransposenabladeltapsi,inc,
		  tmpdeltaresid,inc);
      cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)];
      if(tmpbasis->ecsflag){
	coeff=0.5*Asqval;//*exp(2.*II*ecstheta);
      }
      else{
	coeff=0.5*Asqval;
      }
      cblas_zaxpy(tmporder,&coeff,tmpodeltapsi,inc,
		  tmpdeltaresid,inc);

      //apply tmpdeltaresid to residxt[n]
      cblas_zaxpy(tmporder,&one,tmpdeltaresid,inc,&residxt[n][eltindx],inc);

      delete [] tmpdeltaresid;
    }

    
    delete [] tmpodeltapsi;
    delete [] tmpVdeltapsi;
    delete [] tmpnablasqdeltapsi;
    delete [] tmpnabladeltapsi;
    delete [] tmptransnabla;
    delete [] tmptransposenabladeltapsi;

    

  }
}



////old way for applydeltapsi
//void wf::applydeltapsi(int ncorr,int nmax, cmplx* __restrict deltapsi,
//			 cmplx* __restrict * __restrict guessxt,
//			 cmplx* __restrict * __restrict residxt, 
//			 temporalbasis* tbas ){
//  cmplx II=cmplx(0.,1.);
//  rl ecstheta=gbas->ecstheta;
//  int ot=tbas->order;
//  //apply deltapsi vector to entry n of guessxt, adjust residuals by
//  //L.deltapsi
//
//  cout << "inside applydeltapsi\n";
//  cout <<"umat\n";
//  printmat(tbas->umat,ot,ot);
//  cout <<"qmat\n";
//  printmat(tbas->qmat,ot,ot);
//
//  //apply deltapsi vector to entry n of guessxt
//  //for(int i=0;i<nbf;i++){
//  //  guessxt[ncorr][i]+=deltapsi[i];
//  //}
//  int inc=1;
//  cmplx one=1.;
//  cblas_zaxpy(nbf,&one,deltapsi,inc,guessxt[ncorr],inc);
//
//  //find H0.deltapsi, nabla.deltapsi, Olap.deltapsi
//  //cmplx* __restrict H0deltapsi=gbas->Hpsi(deltapsi);
//  cmplx* __restrict Vdeltapsi=gbas->Vpsi(deltapsi);
//  cmplx* __restrict nablasqdeltapsi=gbas->nablasqpsi(deltapsi);
//  cmplx* __restrict nabladeltapsi=gbas->nablapsi(deltapsi);
//  cmplx* __restrict transposenabladeltapsi=gbas->transposenablapsi(deltapsi);
//  cmplx* __restrict Olapdeltapsi=gbas->Olappsi(deltapsi);
//  cmplx* __restrict nablaolapdeltapsi=gbas->nablapsi(Olapdeltapsi);
//  cmplx* __restrict olapnabladeltapsi=gbas->Olappsi(nabladeltapsi);
//  cmplx* __restrict Olapsqdeltapsi=gbas->Olappsi(Olapdeltapsi);
//
////  //find Amat.Umat
////  cmplx* __restrict AUmat=initializearray(ot*ot);
////  cmplx alpha=1.;
////  cmplx beta=0.;
////  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
////	      tbas->Amat,ot,tbas->umat,ot,&beta,AUmat,ot);
////  cmplx* __restrict UAmat=initializearray(ot*ot);
////  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
////	      tbas->umat,ot,tbas->Amat,ot,&beta,UAmat,ot);
//  
////  //find Asqmat.Umat
////  cmplx* __restrict AsqUmat=initializearray(ot*ot);
//////  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
//////	      tbas->Amat,ot,AUmat,ot,&beta,AsqUmat,ot);
////  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
////	      tbas->Asqmat,ot,tbas->umat,ot,&beta,AsqUmat,ot);
//
//
//  //correct residuals according to L.deltapsi
//  //for tridiagonal temporal basis, umat will have bandwith of at most 2
//  //q will have nonzero elements only for n<ncorr
//  //Amat and Asqmat may have nonzero elements of any order
//  //umat has bandwidth 2
//  //for(int n=max(0,ncorr-2);n<min(nmax,ncorr+2);n++){
//  for(int n=0;n<nmax;n++){
//    cmplx uval=tbas->umat[arrayindx(n,ncorr,ot,ot)];
//    //cout << "applydeltapsi uval "<<ncorr<<" "<<n<<" "<<uval<<"\n";
//    //cmplx coeff=-uval;
//    //cblas_zaxpy(nbf,&coeff,H0deltapsi,inc,residxt[n],inc);
//    cmplx coeff=-uval;
//    cblas_zaxpy(nbf,&coeff,Vdeltapsi,inc,residxt[n],inc);
//  }
//  //qmat elts are nonzero only for n<=ncorr
//  //for(int n=0;n<=max(ncorr,2);n++){
//  for(int n=0;n<nmax;n++){
//    cmplx qval=tbas->qmat[arrayindx(n,ncorr,ot,ot)];
//    //cout << "applydeltapsi qval "<<ncorr<<" "<<n<<" "<<qval<<"\n";
//    cmplx coeff=II*qval;
//    cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
//  }
//  //for(int n=0;n<nmax;n++){
//  //  cmplx AUval=AUmat[arrayindx(n,ncorr,ot,ot)];
//  //  cmplx AsqUval=AsqUmat[arrayindx(n,ncorr,ot,ot)];
//  //  //cmplx coeff=II*AUval;
//  //  cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)];
//  //  cmplx ATval=tbas->Amat[arrayindx(ncorr,n,ot,ot)];//value of A transpose
//  //  cmplx coeff=-II*ATval/2.;
//  //  cblas_zaxpy(nbf,&coeff,nabladeltapsi,inc,residxt[n],inc);
//  //  coeff=II*Aval/2.;
//  //  cblas_zaxpy(nbf,&coeff,transposenabladeltapsi,inc,residxt[n],inc);
//  //  cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)];
//  //  //coeff=-AsqUval;
//  //  coeff=Asqval;
//  //  cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
//  //}
//  for(int n=0;n<nmax;n++){
//    cmplx uval=tbas->umat[arrayindx(n,ncorr,ot,ot)];
//    //cmplx AUval=AUmat[arrayindx(n,ncorr,ot,ot)];
//    //cmplx UAval=UAmat[arrayindx(n,ncorr,ot,ot)];
//    //cmplx AsqUval=AsqUmat[arrayindx(n,ncorr,ot,ot)];
//    //cmplx coeff=II*AUval;
//    cmplx coeff=0.5*uval;
//    cblas_zaxpy(nbf,&coeff,nablasqdeltapsi,inc,residxt[n],inc);
//    cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)];
//    
//    coeff=0.5*II*Aval;
//    cblas_zaxpy(nbf,&coeff,nabladeltapsi,inc,residxt[n],inc);
//    coeff=-0.5*II*Aval;
//    cblas_zaxpy(nbf,&coeff,transposenabladeltapsi,inc,residxt[n],inc);
//
//    cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)];
//    //coeff=-AsqUval;
//    coeff=0.5*Asqval;
//    cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
//  }
//  
//  
//
////   //old way: don't worrry about which elements are nonzero, just add them in
////   //(works slowly because it spends time on adding 0)
////  for(int n=0;n<nmax;n++){
////    cmplx uval=tbas->umat[arrayindx(n,ncorr,ot,ot)];
////    cmplx qval=tbas->qmat[arrayindx(n,ncorr,ot,ot)];
////    cmplx AUval=AUmat[arrayindx(n,ncorr,ot,ot)];
////    cmplx AsqUval=AsqUmat[arrayindx(n,ncorr,ot,ot)];
////    cout << "ncorr, n, uval, qval,AUval,AsqUval\t"<<ncorr<<"\t"<<n<<"\t"<<uval<<"\t"<<qval<<"\t"<<AUval<<"\t"<<AsqUval<<"\n";
////    //for(int i=0;i<nbf;i++){
////    //  residxt[n][i]+=
////    //	(-uval*H0deltapsi[i]//-Umat.(H0.deltapsi)
////    //	+AUval*nabladeltapsi[i]//+Amat.Umat.(nabla.deltapsi)
////    //	-AsqUval*Olapdeltapsi[i]//-1/2 Asqmat.Umat.(Olap.deltapsi)
////    //	 +II*qval*Olapdeltapsi[i]);//+II*qmat.(Olap.deltapsi)
////    //}
////    
////    cmplx coeff=-uval;
////    cblas_zaxpy(nbf,&coeff,H0deltapsi,inc,residxt[n],inc);
////    coeff=II*AUval*0.5;
////    cblas_zaxpy(nbf,&coeff,nablaolapdeltapsi,inc,residxt[n],inc);
////    cblas_zaxpy(nbf,&coeff,olapnabladeltapsi,inc,residxt[n],inc);
////    coeff=(-0.5*AsqUval);
////    cblas_zaxpy(nbf,&coeff,Olapsqdeltapsi,inc,residxt[n],inc);
////    coeff=II*qval;
////    cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
////  }
//  
//  
//
//  //delete [] AUmat;
//  //delete [] AsqUmat;
//  //delete [] H0deltapsi;
//  delete [] Vdeltapsi;
//  delete [] nablasqdeltapsi;
//  delete [] Olapdeltapsi;
//  delete [] Olapsqdeltapsi;
//  delete [] nabladeltapsi;
//  delete [] transposenabladeltapsi;
//  delete [] nablaolapdeltapsi;
//  delete [] olapnabladeltapsi;
//
//}

//apply corrections in the same way that the linear system is set up
void wf::deltapsiresidual_Aterms(int ncorr,int nmax, cmplx* __restrict deltapsi,
				 cmplx Aratio, cmplx Asqratio,
				 cmplx* __restrict * __restrict residxt, 
				 temporalbasis* tbas ){
  cmplx II=cmplx(0.,1.);
  rl ecstheta=gbas->ecstheta;
  int ot=tbas->order;
  //apply deltapsi vector to entry n of guessxt, adjust residuals by
  //L.deltapsi

  //cout << "inside applydeltapsi\n";
  //cout <<"umat\n";
  //printmat(tbas->umat,ot,ot);
  //cout <<"qmat\n";
  //printmat(tbas->qmat,ot,ot);

  //apply deltapsi vector to entry n of guessxt
  //for(int i=0;i<nbf;i++){
  //  guessxt[ncorr][i]+=deltapsi[i];
  //}
  int inc=1;
  cmplx one=1.;
  //cblas_zaxpy(nbf,&one,deltapsi,inc,guessxt[ncorr],inc);


  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
    int tmporder=tmpbasis->order;
    int eltindx=gbas->eltfirstindices[eltnum];
    cmplx* __restrict tmph=tmpbasis->Hmat;
    cmplx* __restrict tmpnablasq=tmpbasis->nablasqmat;
    cmplx* __restrict tmpV=tmpbasis->Vmat;
    cmplx* __restrict tmpo=tmpbasis->omat;
    cmplx* __restrict tmpnabla=tmpbasis->nablamat;
    cmplx* __restrict tmptransnabla=transpose(tmpnabla,tmporder,tmporder);


    //cout << "tmpnabla\n";
    //printmat(tmpnabla,tmporder,tmporder);
    //cout << "tmpnablasq\n";
    //printmat(tmpnablasq,tmporder,tmporder);


    cmplx alpha=1.;
    cmplx beta=0.;

    //find product of Hamiltonian & overlap arrays within temporary basis
    cmplx* __restrict tmpodeltapsi=tmpbasis->Xmatpsi(tmpbasis->omat,
						&deltapsi[eltindx]);
    cmplx* __restrict tmpVdeltapsi=tmpbasis->Xmatpsi(tmpbasis->Vmat,
						&deltapsi[eltindx]);
    cmplx* __restrict tmpnablasqdeltapsi=tmpbasis->Xmatpsi(tmpbasis->nablasqmat,
						&deltapsi[eltindx]);
    cmplx* __restrict tmpnabladeltapsi=tmpbasis->Xmatpsi(tmpbasis->nablamat,
						&deltapsi[eltindx]);

    cmplx* __restrict tmptransposenabladeltapsi=tmpbasis->
      Xmatpsi(tmptransnabla,&deltapsi[eltindx]);
    
    //sum L.psi for current element & apply to residxt
    for(int n=0;n<nmax;n++){
//      cmplx* __restrict tmpdeltaresid=initializearray(tmporder);
//      cmplx coeff=0.;
//      //      if(n<=max(ncorr,2)){
//	cmplx qval=tbas->qmat[arrayindx(n,ncorr,ot,ot)];
//	coeff=II*qval;
//	cblas_zaxpy(tmporder,&coeff,tmpodeltapsi,inc,tmpdeltaresid,inc);
//	//}
//	//if(n>=(max(0,ncorr-2)) and (n<min(nmax,ncorr+2))){
//	cmplx uval=tbas->umat[arrayindx(n,ncorr,ot,ot)];
//	coeff=-uval;
//	cblas_zaxpy(tmporder,&coeff,tmpVdeltapsi,inc,tmpdeltaresid,inc);
//
//	coeff=0.5*uval;
//	cblas_zaxpy(tmporder,&coeff,tmpnablasqdeltapsi,inc,tmpdeltaresid,inc);
//	//}


      //terms involving vector potential
      cmplx coeff=0.;
      cmplx* __restrict tmpdeltaresid=initializearray(tmporder);
      cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)]-
	Aratio*tbas->umat[arrayindx(n,ncorr,ot,ot)];
      if(tmpbasis->ecsflag){
	coeff=0.5*II*Aval;//*exp(II*ecstheta);
      }
      else{
	coeff=0.5*II*Aval;
      }
      //coeff*=-1.;
      cblas_zaxpy(tmporder,&coeff,tmpnabladeltapsi,inc,tmpdeltaresid,inc);

      if(tmpbasis->ecsflag){
	coeff=-0.5*II*Aval;//*exp(-II*ecstheta);//exp(II*ecstheta);
      }
      else{
	coeff=-0.5*II*Aval;
      }
      //coeff*=-1.;
      cblas_zaxpy(tmporder,&coeff,tmptransposenabladeltapsi,inc,
		  tmpdeltaresid,inc);
      cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)]-
	Asqratio*tbas->umat[arrayindx(n,ncorr,ot,ot)];
      if(tmpbasis->ecsflag){
	coeff=0.5*Asqval;//*exp(2.*II*ecstheta);
      }
      else{
	coeff=0.5*Asqval;
      }
      cblas_zaxpy(tmporder,&coeff,tmpodeltapsi,inc,
		  tmpdeltaresid,inc);

      //apply tmpdeltaresid to residxt[n]
      cblas_zaxpy(tmporder,&one,tmpdeltaresid,inc,&residxt[n][eltindx],inc);

      delete [] tmpdeltaresid;
    }

    
    delete [] tmpodeltapsi;
    delete [] tmpVdeltapsi;
    delete [] tmpnablasqdeltapsi;
    delete [] tmpnabladeltapsi;
    delete [] tmptransnabla;
    delete [] tmptransposenabladeltapsi;

    

  }
}


////Old way to find deltapsiresidual_Aterms
//void wf::deltapsiresidual_Aterms(int ncorr,int nmax, cmplx* __restrict deltapsi,
//				 cmplx Aratio, cmplx Asqratio,
//				 cmplx* __restrict * __restrict residxt, 
//				 temporalbasis* tbas ){
//  cmplx II=cmplx(0.,1.);
//  rl ecstheta=gbas->ecstheta;
//  int ot=tbas->order;
//  //apply deltapsi vector to entry n of guessxt, adjust residuals by
//  //L.deltapsi
//
//  //cout << "inside deltapsiresidual_Aterms\n";
//  //cout <<"umat\n";
//  //printmat(tbas->umat,ot,ot);
//  //cout <<"qmat\n";
//  //printmat(tbas->qmat,ot,ot);
//
//  //apply deltapsi vector to entry n of guessxt
//  //for(int i=0;i<nbf;i++){
//  //  guessxt[ncorr][i]+=deltapsi[i];
//  //}
//  int inc=1;
//  cmplx one=1.;
//
//  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
//    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
//    int tmporder=tmpbasis->order;
//    int eltindx=gbas->eltfirstindices[eltnum];
//    cmplx* __restrict tmpo=tmpbasis->omat;
//    cmplx* __restrict tmpnabla=tmpbasis->nablamat;
//    cmplx* __restrict tmptransnabla=transpose(tmpnabla,tmporder,tmporder);
//
//
//    //cout << "tmpnabla\n";
//    //printmat(tmpnabla,tmporder,tmporder);
//    //cout << "tmpnablasq\n";
//    //printmat(tmpnablasq,tmporder,tmporder);
//
//
//    cmplx alpha=1.;
//    cmplx beta=0.;
//
//    //find product of Hamiltonian & overlap arrays within temporary basis
//    cmplx* __restrict tmpodeltapsi=tmpbasis->Xmatpsi(tmpbasis->omat,
//						     &deltapsi[eltindx]);
//    cmplx* __restrict tmposqdeltapsi=tmpbasis->Xmatpsi(tmpbasis->omat,
//						       tmpodeltapsi);
//    cmplx* __restrict tmpnabladeltapsi=tmpbasis->Xmatpsi(tmpbasis->nablamat,
//							 &deltapsi[eltindx]);
//
//    cmplx* __restrict tmptransposenabladeltapsi=tmpbasis->
//      Xmatpsi(tmptransnabla,&deltapsi[eltindx]);
//    
//    //sum L.psi for current element & apply to residxt
//    for(int n=0;n<nmax;n++){
//      cmplx* __restrict tmpdeltaresid=initializearray(tmporder);
//      cmplx coeff=0.;
//
//
//      cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)]-
//	Aratio*tbas->umat[arrayindx(n,ncorr,ot,ot)];
//      if(tmpbasis->ecsflag){
//	coeff=0.5*II*Aval;//*exp(II*ecstheta);
//      }
//      else{
//	coeff=0.5*II*Aval;
//      }
//      //coeff*=-1.;
////      cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)];
////      if(tmpbasis->ecsflag){
////	Aval*=exp(II*ecstheta);
////      }
////      coeff=0.5*II*(Aval-Aratio*tbas->umat[arrayindx(n,ncorr,ot,ot)]);
////      coeff*=-1;
//
//      cblas_zaxpy(tmporder,&coeff,tmpnabladeltapsi,inc,tmpdeltaresid,inc);
//
//      if(tmpbasis->ecsflag){
//	coeff=-0.5*II*Aval;//*exp(-II*ecstheta);
//      }
//      else{
//	coeff=-0.5*II*Aval;
//      }
//      //coeff*=-1.;
//
//
////      coeff=-0.5*II*(Aval-Aratio*tbas->umat[arrayindx(n,ncorr,ot,ot)]);
////      coeff*=-1;
//
//      cblas_zaxpy(tmporder,&coeff,tmptransposenabladeltapsi,inc,
//		  tmpdeltaresid,inc);
//      cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)]-
//	Asqratio*tbas->umat[arrayindx(n,ncorr,ot,ot)];
//      if(tmpbasis->ecsflag){
//	coeff=0.5*Asqval;//*exp(2.*II*ecstheta);
//      }
//      else{
//	coeff=0.5*Asqval;
//      }
//      coeff*=-1.;//*pow(2.,-.5);
//
//
////      cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)];
////      if(tmpbasis->ecsflag){
////	Asqval*=exp(2.*II*ecstheta);
////      }
////      coeff=0.5*(Asqval-Asqratio*tbas->umat[arrayindx(n,ncorr,ot,ot)]);
////      coeff*=-1.;
//
//      cblas_zaxpy(tmporder,&coeff,tmpodeltapsi,inc,
//		  tmpdeltaresid,inc);
//
//      //apply tmpdeltaresid to residxt[n]
//      cblas_zaxpy(tmporder,&one,tmpdeltaresid,inc,&residxt[n][eltindx],inc);
//
//      delete [] tmpdeltaresid;
//    }
//
//    
//    delete [] tmpodeltapsi;
//    delete [] tmposqdeltapsi;
//    delete [] tmpnabladeltapsi;
//    delete [] tmptransnabla;
//    delete [] tmptransposenabladeltapsi;
//
//    
//
//  }
//}




//////////////////////////////////////////////////////////////////////////







//void wf::deltapsiresidual_Aterms(int ncorr,int nmax, cmplx* __restrict deltapsi,
//				 cmplx Aratio, cmplx Asqratio,
//				 cmplx* __restrict * __restrict guessxt,
//				 cmplx* __restrict * __restrict residxt, 
//				 temporalbasis* tbas ){
//  //finds residual resulting from those terms of the Hamiltonian which weren't
//  //included in deltawf_eigenvector -- namely, (Amat-Aval*Umat) and
//  //(Asqmat-Asqval*Umat)
//
//  cmplx II=cmplx(0.,1.);
//  int ot=tbas->order;
//  //apply deltapsi vector to entry n of guessxt, adjust residuals by
//  //L.deltapsi
//
//  //cout << "inside applydeltapsi\n";
//  //cout <<"umat\n";
//  //printmat(tbas->umat,ot,ot);
//  //cout <<"qmat\n";
//  //printmat(tbas->qmat,ot,ot);
//
//  //apply deltapsi vector to entry n of guessxt
//  //for(int i=0;i<nbf;i++){
//  //  guessxt[ncorr][i]+=deltapsi[i];
//  //}
//  int inc=1;
//  cmplx one=1.;
//  cblas_zaxpy(nbf,&one,deltapsi,inc,guessxt[ncorr],inc);
//
//  //find H0.deltapsi, nabla.deltapsi, Olap.deltapsi
//  //cmplx* __restrict H0deltapsi=gbas->Hpsi(deltapsi);
//  cmplx* __restrict nabladeltapsi=gbas->nablapsi(deltapsi);
//  //cmplx* __restrict transposenabladeltapsi=gbas->transposenablapsi(deltapsi);
//  cmplx* __restrict Olapdeltapsi=gbas->Olappsi(deltapsi);
//  //cmplx* __restrict nablaolapdeltapsi=gbas->nablapsi(Olapdeltapsi);
//  //cmplx* __restrict olapnabladeltapsi=gbas->Olappsi(nabladeltapsi);
//  //cmplx* __restrict Olapsqdeltapsi=gbas->Olappsi(Olapdeltapsi);
//
////  //find Amat.Umat
////  cmplx* __restrict AUmat=initializearray(ot*ot);
////  cmplx alpha=1.;
////  cmplx beta=0.;
////  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
////	      tbas->Amat,ot,tbas->umat,ot,&beta,AUmat,ot);
////  //find Asqmat.Umat
////  cmplx* __restrict AsqUmat=initializearray(ot*ot);
//////  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
//////	      tbas->Amat,ot,AUmat,ot,&beta,AsqUmat,ot);
////  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
////	      tbas->Asqmat,ot,tbas->umat,ot,&beta,AsqUmat,ot);
//
//
//  //correct residuals according to L.deltapsi
//  //for tridiagonal temporal basis, umat will have bandwith of at most 2
//  //q will have nonzero elements only for n<ncorr
//  //Amat and Asqmat may have nonzero elements of any order
//  //umat has bandwidth 2
////  for(int n=max(0,ncorr-2);n<min(nmax,ncorr+2);n++){
////    cmplx uval=tbas->umat[arrayindx(n,ncorr,ot,ot)];
////    //cout << "applydeltapsi uval "<<ncorr<<" "<<n<<" "<<uval<<"\n";
////    cmplx coeff=-uval;
////    cblas_zaxpy(nbf,&coeff,H0deltapsi,inc,residxt[n],inc);
////  }
////  //qmat elts are nonzero only for n<=ncorr
////  for(int n=0;n<=max(ncorr,1);n++){
////    cmplx qval=tbas->qmat[arrayindx(n,ncorr,ot,ot)];
////    //cout << "applydeltapsi qval "<<ncorr<<" "<<n<<" "<<qval<<"\n";
////    cmplx coeff=II*qval;
////    cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
////  }
//  //for(int n=0;n<nmax;n++){
//  //  cmplx AUval=AUmat[arrayindx(n,ncorr,ot,ot)];
//  //  cmplx AsqUval=AsqUmat[arrayindx(n,ncorr,ot,ot)];
//  //  //cmplx coeff=II*AUval;
//  //  cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)];
//  //  cmplx ATval=tbas->Amat[arrayindx(ncorr,n,ot,ot)];//value of A transpose
//  //  cmplx coeff=-II*ATval/2.;
//  //  cblas_zaxpy(nbf,&coeff,nabladeltapsi,inc,residxt[n],inc);
//  //  coeff=II*Aval/2.;
//  //  cblas_zaxpy(nbf,&coeff,transposenabladeltapsi,inc,residxt[n],inc);
//  //  cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)];
//  //  //coeff=-AsqUval;
//  //  coeff=Asqval;
//  //  cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
//  //}
//  for(int n=0;n<nmax;n++){
//    //cmplx AUval=AUmat[arrayindx(n,ncorr,ot,ot)];
//    //cmplx AsqUval=AsqUmat[arrayindx(n,ncorr,ot,ot)];
//    //cmplx coeff=II*AUval;
//    cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)]-
//      Aratio*triD_U(n,ncorr);
//    cmplx coeff=-II*Aval;
//    cblas_zaxpy(nbf,&coeff,nabladeltapsi,inc,residxt[n],inc);
//    cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)]-
//      Asqratio*triD_U(n,ncorr);
//    //coeff=-AsqUval;
//    coeff=0.5*Asqval;
//    cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
//  }
//  
//  
//
////   //old way: don't worrry about which elements are nonzero, just add them in
////   //(works slowly because it spends time on adding 0)
////  for(int n=0;n<nmax;n++){
////    cmplx uval=tbas->umat[arrayindx(n,ncorr,ot,ot)];
////    cmplx qval=tbas->qmat[arrayindx(n,ncorr,ot,ot)];
////    cmplx AUval=AUmat[arrayindx(n,ncorr,ot,ot)];
////    cmplx AsqUval=AsqUmat[arrayindx(n,ncorr,ot,ot)];
////    cout << "ncorr, n, uval, qval,AUval,AsqUval\t"<<ncorr<<"\t"<<n<<"\t"<<uval<<"\t"<<qval<<"\t"<<AUval<<"\t"<<AsqUval<<"\n";
////    //for(int i=0;i<nbf;i++){
////    //  residxt[n][i]+=
////    //	(-uval*H0deltapsi[i]//-Umat.(H0.deltapsi)
////    //	+AUval*nabladeltapsi[i]//+Amat.Umat.(nabla.deltapsi)
////    //	-AsqUval*Olapdeltapsi[i]//-1/2 Asqmat.Umat.(Olap.deltapsi)
////    //	 +II*qval*Olapdeltapsi[i]);//+II*qmat.(Olap.deltapsi)
////    //}
////    
////    cmplx coeff=-uval;
////    cblas_zaxpy(nbf,&coeff,H0deltapsi,inc,residxt[n],inc);
////    coeff=II*AUval*0.5;
////    cblas_zaxpy(nbf,&coeff,nablaolapdeltapsi,inc,residxt[n],inc);
////    cblas_zaxpy(nbf,&coeff,olapnabladeltapsi,inc,residxt[n],inc);
////    coeff=(-0.5*AsqUval);
////    cblas_zaxpy(nbf,&coeff,Olapsqdeltapsi,inc,residxt[n],inc);
////    coeff=II*qval;
////    cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
////  }
//  
//  
//
//  //delete [] AUmat;
//  //delete [] AsqUmat;
//  //delete [] H0deltapsi;
//  delete [] Olapdeltapsi;
//  //delete [] Olapsqdeltapsi;
//  delete [] nabladeltapsi;
//  //delete [] transposenabladeltapsi;
//  //delete [] nablaolapdeltapsi;
//  //delete [] olapnabladeltapsi;
//
//}







cmplx* __restrict 
wf::temporalvectorprojection_revec(int ot, cmplx* __restrict levec,
			     cmplx* __restrict * __restrict psixt){
  //given psixt and a left eigenvector levec, find the projection of psixt
  //onto levec's corresponding right eigenvector
  
  int nt=ot-1;

  cmplx* __restrict U_red=initializearray((ot-1)*(ot-1));
  for(int i=1;i<ot;i++){
    for(int j=1;j<ot;j++){
      U_red[arrayindx(i-1,j-1,ot-1,ot-1)]=triD_U(i,j);
    }
  }
  //find (vec)**H.U_red = (U_red**H.vec)**H
  cmplx alpha=1.;
  cmplx beta=0.;
  int incx=1;
  int incy=1;
  cmplx* __restrict Ulevec=initializearray(ot);
  cblas_zgemv(CblasColMajor,CblasConjTrans,nt,nt,&alpha,U_red,
	      nt,levec+1,incx,&beta,Ulevec+1,incy);
  for(int n=0;n<ot;n++){
    Ulevec[n]=conj(Ulevec[n]);
  }
  //cout << "levec\n";
  //printmat(levec,ot,1);
  //cout << "Ulevec\n";
  //printmat_scipy(Ulevec,ot,1);

  cmplx* __restrict retvec=initializearray(nbf);
  //use i as inner loop for easier vectorization
  for(int n=1;n<ot;n++){
    for(int i=0;i<nbf;i++){
      retvec[i]+=Ulevec[n]*psixt[n][i];
    }
  }

  delete [] U_red;
  delete [] Ulevec;
  return retvec;
}








cmplx* __restrict 
wf::temporalvectorprojection_levec(int ot, cmplx* __restrict revec,
			     cmplx* __restrict * __restrict psixt){
  //given psixt and a right eigenvector revec, find the projection of psixt
  //onto levec's corresponding right eigenvector
  
  int nt=ot-1;

  cmplx* __restrict U_red=initializearray((ot-1)*(ot-1));
  for(int i=1;i<ot;i++){
    for(int j=1;j<ot;j++){
      U_red[arrayindx(i-1,j-1,ot-1,ot-1)]=triD_U(i,j);
    }
  }
  //find U_red.revec
  cmplx alpha=1.;
  cmplx beta=0.;
  int incx=1;
  int incy=1;
  cmplx* __restrict Urevec=initializearray(ot);
  cblas_zgemv(CblasColMajor,CblasNoTrans,nt,nt,&alpha,U_red,
	      nt,revec+1,incx,&beta,Urevec+1,incy);
  cout << "revec\n";
  printmat(revec,ot,1);
  cout << "Urevec\n";
  printmat_scipy(Urevec,ot,1);

  cmplx* __restrict retvec=initializearray(nbf);
  //use i as inner loop for easier vectorization
  for(int n=1;n<ot;n++){
    for(int i=0;i<nbf;i++){
      retvec[i]+=conj(Urevec[n])*psixt[n][i];
    }
  }

  delete [] U_red;
  delete [] Urevec;
  return retvec;
}





cmplx* __restrict 
wf::temporalvector_rhstarget(int ot, cmplx* __restrict levec,
			     cmplx* __restrict * __restrict psixt){
  //given psixt and a left eigenvector levec, find the projection of psixt
  //which should be the rhs in the linear solve.
  //left and right eigenvectors are normalized s.t. W*_n U_red v_m = delta_nm, 
  //where W is left eigenvector and V is right eigenvector.
  //then, linear system will solve for
  //(since Q v_n=lambda U v_n
  //(i O lambda_n U - H U)v_n dpsi_n=R
  //(i O lambda_n - H) v_n dpsi_n = Uinv.R
  //w*_m U (i O lambda_n - H) v_n dpsi_n=w*_m.U.Uinv.R
  //thus, rhs should be w*_m.R
  
  

  cmplx* __restrict retvec=initializearray(nbf);
  //use i as inner loop for easier vectorization
  int inc=1;
  for(int n=1;n<ot;n++){
    cmplx lval=conj(levec[n]);
    cblas_zaxpy(nbf,&lval,psixt[n],inc,retvec,inc);
    //for(int i=0;i<nbf;i++){
    //  retvec[i]+=lval*psixt[n][i];
    //}
  }

  return retvec;
}














void wf::applycorrection_temporalevec(int ot,cmplx* __restrict dpsin, 
				      cmplx* __restrict vec,
				      cmplx* __restrict * __restrict guessxt,
				      cmplx* __restrict * __restrict residxt,
				      temporalbasis* tbas ){

  for(int n=0;n<ot;n++){
    cmplx* __restrict tmpvec=arraymultiply(dpsin,nbf,vec[n]);
    applydeltapsi(n,ot,tmpvec,guessxt,residxt,tbas);
    delete [] tmpvec;
  }

}

void wf::initialguess(int nmax, cmplx* __restrict psi0, 
		      cmplx* __restrict * __restrict psixt,
		      cmplx* __restrict * __restrict residxt,
		      temporalbasis* tbas){
  //choose psixt s.t. psixt[0]=psi0, Sum_m Unm.psixt[m]=0 for n>0.  This will
  //help with numerical stability, since it will mean that H.U.psixt=0 for
  //n>0, so that there's no repeated multiplication by Hamiltonian matrix to
  //worry about.

  //cmplx* vec=ICvec_(nmax);
  rl dt=tbas->tmax-tbas->tmin;
  for(int n=0;n<nmax;n++){
    zeroarray(psixt[n],nbf);
    zeroarray(residxt[n],nbf);
  }
  //set level 0 according to initial condition
  applydeltapsi(0,nmax,psi0,psixt,residxt,tbas);
  //zeroarray(residxt[0],nbf);
  //set initial guess to be a constant function
  applydeltapsi(1,nmax,psi0,psixt,residxt,tbas);
  //set level 1 so as to yield zero residual on level 0
  //cmplx* __restrict psi1=residualcorrection(1,0,residxt[0],tbas);
  //applydeltapsi(1,nmax,psi1,psixt,residxt,tbas);
  //delete [] psi1;
}

cmplx * __restrict * __restrict 
wf::deltawfxt_temporaleigenvectors(int nmax,int recursionlvl, 
				   int maxrecursionlvl, rl accuracygoal,
				   cmplx& stepaction,
				   cmplx* __restrict psi0, 
				   cmplx* __restrict * __restrict targetxt, 
				   temporalbasis* tbas){
  //given psi0, targetxt: 
  //1)impose initial conditions & make 1st level
  //correction so as to have zero residual on level 0

  //2)project targetxt-residxt onto temporal eigenvectors 

  //3)for each temporal eigenvector, calculate corresponding correction to
  //psixt

   //4)if deviation between residual and target is too large, call self
  //recursively, with zero initial condition
  
  int ot=tbas->order;


  cmplx* __restrict * __restrict guessxt=new cmplx*[nmax];
  cmplx* __restrict * __restrict residxt=new cmplx*[nmax];
  //cmplx* __restrict * __restrict deltapsixt=new cmplx*[nmax];
  //cmplx* __restrict * __restrict deltadiffxt=new cmplx*[nmax];
  //  cmplx* __restrict * __restrict dpsi=new cmplx*[nmax];
  //cmplx* __restrict * __restrict diff_rvecproj=new cmplx*[nmax];
  //cmplx* __restrict * __restrict nexttarget=new cmplx*[nmax];
  for(int n=0;n<nmax;n++){
    guessxt[n]=initializearray(nbf);
    residxt[n]=initializearray(nbf);
    //deltapsixt[n]=initializearray(nbf);
    //deltadiffxt[n]=initializearray(nbf);
    //nexttarget[n]=initializearray(nbf);
  }
  //set up initial guess: diffxt will be action difference between guess & 
  //target
  //if(recursionlvl==0){

  //set up Q,U,UQ and UU matrices for use in generalized eigenvector
  //decomposition
  cmplx* __restrict tmpU=triD_Umat(ot);
  cmplx* __restrict tmpQ=triD_Qmat(ot);
  cmplx* __restrict tmpUU=initializearray(ot*ot);
  cmplx* __restrict tmpUQ=initializearray(ot*ot);

//  //set up matrix to swap zeroth and last temporal coeffs
//  cmplx* __restrict swapmat=initializearray(ot*ot);
//  for(int n=1;n<ot-1;n++){
//    swapmat[arrayindx(n,n,ot,ot)]=1.;
//  }
//  swapmat[arrayindx(0,ot-1,ot,ot)]=1.;
//  swapmat[arrayindx(ot-1,0,ot,ot)]=1.;
//  cmplx* __restrict tmpSU=initializearray(ot*ot);
//  cmplx* __restrict tmpSQ=initializearray(ot*ot);
//
//  cmplx alpha=1.;
//  cmplx beta=0.;
//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
//	      tmpU,ot,tmpQ,ot,&beta,tmpUQ,ot);
//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
//	      tmpU,ot,tmpU,ot,&beta,tmpUU,ot);
//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
//	      swapmat,ot,tmpQ,ot,&beta,tmpSQ,ot);
//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
//	      swapmat,ot,tmpU,ot,&beta,tmpSU,ot);
  
  //cout << "swapmat\n";
  //printmat(swapmat,ot,ot);
  //cout << "tmpSU\n";
  //printmat(tmpSU,ot,ot);

  //for(int n=0;n<ot;n++){
  //  cout << "targetxt\t"<<n<<"\t"<<nthcomponentsignificance(targetxt,n)<<"\n";
  //}


  initialguess(nmax,psi0,guessxt,residxt,tbas);
  cmplx* __restrict * __restrict diffxt=psixtdiff(targetxt,residxt,nmax);
//  if(recursionlvl==0){
//    //set level 1 so as to yield zero residual on level 0
//    cmplx* __restrict psi1=residualcorrection(1,0,diffxt[0],tbas);
//    applydeltapsi(1,nmax,psi1,guessxt,residxt,tbas);
//    delete [] psi1;
//  }
  //for(int n=0;n<ot;n++){
  //  cout << "diffxt\t"<<n<<"\t"<<nthcomponentsignificance(diffxt,n)<<"\n";
  //}

  cmplx* __restrict * __restrict Udiffxt=Upsi(ot,diffxt);
  //cmplx* __restrict * __restrict Sdiffxt=psi_swaporders(ot,0,ot-1,diffxt);
  //for(int n=0;n<ot;n++){
  //  cout << "Udiffxt\t"<<n<<"\t"<<nthcomponentsignificance(Udiffxt,n)<<"\n";
  //}
  //for(int n=0;n<ot;n++){
  //  cout << "Sdiffxt\t"<<n<<"\t"<<nthcomponentsignificance(Sdiffxt,n)<<"\n";
  //}


  //cancel Udiffxt residuals for n>=1
//  cancelresiduals_generalized_eigenvector(
//	nbf,ot,tmpUQ,tmpUU,guessxt,residxt,Udiffxt,tbas);

//  //cancel Udiffxt residuals for n>=1 using Olap preconditioner
//  cancelresiduals_generalized_eigenvector_Olap_preconditioner(
//	nbf,ot,tmpUQ,tmpUU,guessxt,residxt,Udiffxt,tbas);

//  cancelresiduals_generalized_eigenvector_Olap_preconditioner(nbf,ot,tmpSQ,tmpSU,guessxt,residxt,
//					  Sdiffxt,tbas);
  //cout << "canceled Udiffxt\n\n";
  //for(int n=0;n<ot;n++){
  //  cout << "residxtxt\t"<<n<<"\t"<<nthcomponentsignificance(residxt,n)<<"\n";
  //}
  
//  //delete old diffxt, Udiffxt
//  for(int n=0;n<ot;n++){
//    delete [] diffxt[n];
//    delete [] Udiffxt[n];
//  }
//  delete [] diffxt;
//  delete [] Udiffxt;
//
//  //find new diffxt
//  diffxt=psixtdiff(targetxt,residxt,nmax);
//  //diffxt=psixtdiff(residxt,targetxt,nmax);
//  Udiffxt=Upsi(ot,diffxt);
  //for(int n=0;n<ot;n++){
  //  cout << "new diffxt\t"<<n<<"\t"<<nthcomponentsignificance(diffxt,n)<<"\n";
  //}
  //for(int n=0;n<ot;n++){
  //  cout << "new Udiffxt\t"<<n<<"\t"<<nthcomponentsignificance(Udiffxt,n)<<"\n";
  //}

  //cancel diffxt residuals for n>=1
  cancelresiduals_generalized_eigenvector(nbf,ot,tmpQ,tmpU,guessxt,residxt,
					  diffxt,tbas);


//  //find the residual remaining due to the fact that Amat!=Aratio*Umat
//  cmplx Aratio=tbas->Amat[arrayindx(1,1,ot,ot)]/
//    tbas->umat[arrayindx(1,1,ot,ot)];
//  cmplx Asqratio=tbas->Asqmat[arrayindx(1,1,ot,ot)]/
//    tbas->umat[arrayindx(1,1,ot,ot)];
//  for(int n=1;n<nmax;n++){
//    deltapsiresidual_Aterms(n,nmax,guessxt[n],Aratio,Asqratio,nexttarget,tbas);
//  }
  
  //convert guessxt to legendre basis
  cmplx* __restrict * __restrict guessxt_legendre=
    triD_legendre_conversion(nmax,guessxt);

  



  //delete old diffxt, Udiffxt
  for(int n=0;n<ot;n++){
    delete [] diffxt[n];
    delete [] Udiffxt[n];
  }
  delete [] diffxt;
  delete [] Udiffxt;

  //find new diffxt
  diffxt=psixtdiff(targetxt,residxt,nmax);
  //diffxt=psixtdiff(residxt,targetxt,nmax);
  Udiffxt=Upsi(ot,diffxt);
  //for(int n=0;n<ot;n++){
  //  cout << "even newer diffxt\t"<<n<<"\t"<<nthcomponentsignificance(diffxt,n)<<"\n";
  //}
  //for(int n=0;n<ot;n++){
  //  cout << "even newer Udiffxt\t"<<n<<"\t"<<nthcomponentsignificance(Udiffxt,n)<<"\n";
  //}

//    //What corrections are needed at level 0 to cancel the remaining residuals?
//  cmplx* __restrict dpsi0=residualcorrection(0,0,diffxt[0],tbas);
//  cmplx* __restrict dpsi0minuspsi0=arraydiff(-1.,dpsi0,1.,psi0,nbf);
//  cout << "dpsi0 significance\t"<<sqrt(abs(gbas->norm(dpsi0)))<<"\n";
//  cout << "dpsi0minuspsi0 significance\t"<<sqrt(abs(gbas->norm(dpsi0minuspsi0)))<<"\n";
  //applydeltapsi(0,ot,dpsi0minuspsi0,guessxt,residxt,tbas);


  //for(int n=0;n<ot;n++){
  //  cout << "residxt\t"<<n<<"\t"<<nthcomponentsignificance(residxt,n)<<"\n";
  //}
  
  //nexttarget is difference between residxt and targetxt
  cmplx* __restrict * __restrict nexttarget=psixtdiff(targetxt,residxt,nmax);
  //cmplx* __restrict * __restrict nexttarget=psixtdiff(residxt,targetxt,nmax);
  cmplx* __restrict * __restrict targetdiff=psixtdiff(nexttarget,targetxt,nmax);
  //for(int n=0;n<ot;n++){
  //  cout << "guessxt\t"<<n<<"\t"<<nthcomponentsignificance(guessxt,n)<<"\n";
  //}
  //for(int n=0;n<ot;n++){
  //  cout << "nexttarget\t"<<n<<"\t"<<nthcomponentsignificance(nexttarget,n)<<"\n";
  //}

  rl dt=tbas->tmax-tbas->tmin;
  rl normsum=normsum_legendre(nmax,dt,guessxt_legendre);
  cout << "recursionlevel, normsum(guessxt)\t"<<recursionlvl<<"\t"<<
      normsum<<"\n";

  cout << "recursionlevel, maxdelta(guessxt)\t"<<recursionlvl<<"\t"<<
      maxdelta(guessxt,ot)<<"\n";
  cout << "recursionlevel, maxdelta(residxt)\t"<<recursionlvl<<"\t"<<
      maxdelta(residxt,ot)<<"\n";
  cout << "recursionlevel, maxdelta(diffxt)\t"<<recursionlvl<<"\t"<<
      maxdelta(diffxt,ot)<<"\n";
  cout << "recursionlevel, maxdelta(nexttarget)\t"<<recursionlvl<<"\t"<<
      maxdelta(nexttarget,ot)<<"\n";
  cout << "recursionlevel, maxdelta(targetdiff)\t"<<recursionlvl<<"\t"<<
      maxdelta(targetdiff,ot)<<"\n";
  cmplx* __restrict lvlpsif=psif(guessxt,ot);
  rl lvlpsifsignificance=sqrt(abs(gbas->norm(lvlpsif)));
  cout << "recursionlevel, lvlpsif significance\t"<<lvlpsifsignificance<<"\n";
  
  cmplx correctionaction=residualaction(nbf,ot,guessxt,nexttarget);
  cout << "recursionlevel, correction action\t"<<abs(correctionaction)<<"\n";

  //if(maxdelta(targetdiff,ot)>pow(accuracygoal,1.5)){
  //if(abs(correctionaction)>pow(accuracygoal,2.)){

  if(normsum>pow(accuracygoal,2)){
  //if(lvlpsifsignificance>pow(accuracygoal,2) and (recursionlvl<10)){
  //if(recursionlvl<5){
    cmplx* __restrict nextpsi0=initializearray(nbf);//dpsi0;//dpsi0minuspsi0;//dpsi0;//initializearray(nbf);
    cmplx* __restrict * __restrict nextdeltaxt=
	deltawfxt_temporaleigenvectors(nmax,recursionlvl+1,maxrecursionlvl,
				       accuracygoal,stepaction,
				       nextpsi0,nexttarget,tbas);
      for(int n=0;n<ot;n++){
	applydeltapsi(n,ot,nextdeltaxt[n],guessxt,residxt,tbas);
      }
      
      //delete temporary vectors
      for(int n=0;n<ot;n++){
	delete [] nextdeltaxt[n];
      }
      delete [] nextdeltaxt;
      delete [] nextpsi0;
    }


  //cout << "calling targetdiff\n";
  cmplx* __restrict * __restrict finalresid=psixtdiff(targetxt,residxt,nmax);
  //zeroarray(finalresid[0],nbf);
  //for(int n=0;n<ot;n++){
  //  cout << "finalresid\t"<<n<<"\t"<<nthcomponentsignificance(finalresid,n)<<"\n";
  //}
  //cout << "calling action\n";
  stepaction=residualaction(nbf,ot,guessxt,finalresid);
  cout << "action due to imperfect solution\t"<<stepaction<<"\n";
    


    //delete temporary arrays
    for(int n=0;n<ot;n++){
      delete [] residxt[n];
      delete [] diffxt[n];
      delete [] Udiffxt[n];
      //delete [] Sdiffxt[n];
      delete [] nexttarget[n];
      delete [] targetdiff[n];
      delete [] finalresid[n];
      delete [] guessxt_legendre[n];
      //delete [] deltadiffxt[n];
    }


    //delete [] deltadiffxt;
    delete [] lvlpsif;
    delete [] residxt;
    delete [] diffxt;
    delete [] Udiffxt;
    delete [] guessxt_legendre;
    //delete [] dpsi0;
    //delete [] dpsi0minuspsi0;
    //delete [] Sdiffxt;
    delete [] nexttarget;
    delete [] targetdiff;
    delete [] finalresid;

    delete [] tmpU;
    delete [] tmpQ;
    delete [] tmpUU;
    delete [] tmpUQ;
    //delete [] swapmat;
    //delete [] tmpSU;
    //delete [] tmpSQ;

    
    

//    for(int n=0;n<ot;n++){
//      arrayscale(guessxt[n],nbf,-1.);
//    }

    return guessxt;
}

void wf::timestep_temporaleigenvectors(rl tmin, rl tmax, potential* Pot,
				       pulse* Pls,wf* sourcewf, rl sourceEn,
				       rl accuracygoal, rl& newdt){
  //set up tridiagonal temporal basis (see triDbasis.cpp)
  rl dt=tmax-tmin;
  temporalbasis* tbas=gbas->tbasis;
  tbas->updateTmats(tmin,tmax,Pls);
  sourcewf->gbas->tbasis->updateTmats(tmin,tmax,Pls);
  cmplx* tbas_triD_basischange=triD_basischange(tbas->order);
  temporalbasis* tbas_triD=new temporalbasis(tbas,
                                         tbas_triD_basischange,Pls,tmin,tmax);
  tbas_triD->updateTmats(tmin,tmax,Pls);

  int ot=tbas->order;//order of temporal basis


  //cout << "tbas_triD matrices\n";
  //cout << "umat\n";
  //printmat_scipy(tbas_triD->umat,ot,ot);
  //cout << "qmat\n";
  //printmat_scipy(tbas_triD->qmat,ot,ot);
  //cout << "Amat\n";
  //printmat_scipy(tbas_triD->Amat,ot,ot);
  //cout << "Asqmat\n";
  //printmat_scipy(tbas_triD->Asqmat,ot,ot);


  //targetxt is the action due to (delta H).sourcewf
  cmplx* __restrict * __restrict targetxt=sourcewf->
    inhomogeneous_action(sourceEn,tbas_triD);
  //cmplx* __restrict * __restrict targetxt=sourcewf->
  //  inhomogeneous_action_lengthgauge(sourceEn,Pls,tbas_triD);
  

  cmplx stepaction=0.;
  cmplx* __restrict * __restrict guessxt=
    deltawfxt_temporaleigenvectors(ot,0,ot,accuracygoal,stepaction,psi,targetxt,tbas_triD);

  //convert guessxt to legendre basis
  cmplx* __restrict * __restrict guessxt_legendre=
    triD_legendre_conversion(ot,guessxt);

  for(int nt=0;nt<ot;nt++){
    cout << nt<<"\tcomponent significance\t"<<
      nthcomponentsignificance(guessxt,nt)<<"\n";
  }

  for(int nt=0;nt<ot;nt++){
    cout << nt<<"\tguessxt_legendre component significance\t"<<
      nthcomponentsignificance(guessxt_legendre,nt)<<"\n";
  }


  //rl truncationerror=lastcomponentsignificance(guessxt,ot);
  //rl truncationerror=nthcomponentsignificance(guessxt,ot-1);
  rl truncationerror=lastlegendrecomponentsignificance(ot,dt,guessxt_legendre);

  cout <<"truncation error\t"<<dt<<"\t"<<
    truncationerror<<"\n";
  rl actionaccumulation=abs(stepaction);
  cout << "action accumulation\t"<<actionaccumulation<<"\n";
  //try new time step=dt*(accuracygoal/truncationerror)**(1./(bfordergoal-1))
  rl error=abs(truncationerror);//actionaccumulation;//truncationerror;//max(truncationerror,maxresidual);
  newdt=dt*pow(abs(accuracygoal/error),(1./(ot-1.)));
  cout << "newdt\t"<<dt<<"\t"<<abs(accuracygoal/error)<<"\t"<<pow(abs(accuracygoal/error),(1./(ot-1.)))<<"\t"<<newdt<<"\n";

  
  //maximum timestep is determined by error in source term
  rl maxnewdt=pow(accuracygoal*factorial(ot)/pow(abs(sourceEn),ot),
		  1./real(ot));
  cout << "maxnewdt\t"<<maxnewdt<<"\n";
  newdt=min(newdt,dt*2.);
  newdt=min(newdt,maxnewdt);
  cout << "newdt\t"<<newdt<<"\n";

  delete [] psi;
  cout << "deleted psi\n";
  psi=psif(guessxt,ot);

  cout << "found new psi\n";

  delete [] tbas_triD_basischange;
  delete tbas_triD;
  cout << "deleted tbas_triD\n";
  for(int n=0;n<ot;n++){
    //delete [] Uguessxt[n];
    delete [] guessxt[n];
    //delete [] residxt[n];
    delete [] targetxt[n];
    delete [] guessxt_legendre[n];
  }
  //delete [] Uguessxt;
  delete [] guessxt;
  //delete [] residxt;
  delete [] targetxt;
  delete [] guessxt_legendre;
  cout << "returning from timestep_temporaleigenvectors\n";
}