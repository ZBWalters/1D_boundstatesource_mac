#include "./classes.h"
#include "./linearsystemsetup.h"
//#include "./temporalbasis.h"
//set up the linear system for the least action propagator



cmplx** tmatsetup_ics(int ot,cmplx* tmat, cmplx* umat,cmplx* icbfmat){
  //icbfmat is basis function matrix corresponding to initial condition basis
  //functions (ie, leftborderfunctionbasis).  Here the convention is that the
  //zeroth basis function has value 1 at the initial time, and is orthogonal
  //to all other basis functions, all of which have value 0 at the initial
  //time.
  
  //Specifying an initial condition is thus equivalent to choosing the
  //coefficients of the zeroth temporal basis function, while solving the
  //least action linear system involves solving for the (ot-1) coefficients
  //(per spatial basis function) which are not specified by the initial
  //condition.  So there's a tricky step which must be performed in order to
  //make the linear system solvable: we must project the residuals generated
  //by the initial condition from a space with ot unknowns per spatial basis
  //function to one with (ot-1) unknowns.  This will be accomplished by
  //projecting the residuals onto the (ot-1) basis functions which are zero at
  //the initial time.  So that the solutions to the linear system are
  //unchanged, we will do this to both sides of the linear system we're
  //setting up.  The net result will be the same as that which would be
  //accomplished using Lagrange multipliers to specify initial conditions,
  //only with fewer unknowns to be solved for.

  //This subroutine returns the temporal matrices necessary to set up such a
  //linear system. It returns one (1 x (ot-1)) matrix to map the initial
  //condition into the subspace of coefficients to be solved for, and one
  //((ot-1) x (ot-1)) matrix which maps the coefficients to be solved for into
  //the same space.  If IC is the (1 x ot) matrix corresponding to the initial
  //condition border function and BF is the ((ot-1) x ot) matrix corresponding
  //to the other basis functions, and U is the (ot x ot) temporal overlap
  //matrix, then 
  //retmat0=IC.U.(BF)^T, and
  //retmat1=BF.U.(BF)^T

  cmplx* ICmat=initializearray(ot);
  int vecnum=0;
  for(int i=0;i<ot;i++){
    ICmat[i]=icbfmat[arrayindx(vecnum,i,ot,ot)];
  }

  cmplx* BFmat=initializearray((ot-1)*ot);
  for(int vecnum=1;vecnum<ot;vecnum++){
    for(int i=0;i<ot;i++){
      BFmat[arrayindx(vecnum-1,i,ot-1,ot)]=icbfmat[arrayindx(vecnum,i,ot,ot)];
    }
  }

  //retmat0=IC.U.(BF)^T
  //calculate using BF.U^(T).IC^(T), so that multiplication only involves 2
  //calls to zgemv
  cmplx alpha=1.;
  cmplx beta=0.;
  cmplx* ICtmpmat=new cmplx[ot];
  cmplx* retmat0=new cmplx[ot-1];
  cblas_zgemv(CblasColMajor,CblasTrans,ot,ot,&alpha,umat,ot,ICmat,1,
	      &beta,ICtmpmat,1);
  cblas_zgemv(CblasColMajor,CblasNoTrans,ot-1,ot,&alpha,BFmat,ot-1,ICtmpmat,1,
	      &beta,retmat0,1);
  delete [] ICtmpmat;
  delete [] ICmat;

  //retmat1=BF.U.(BF)^T
  cmplx* BFtmpmat=new cmplx[(ot-1)*ot];
  cmplx* retmat1=new cmplx[(ot-1)*(ot-1)];
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot-1,ot,ot,&alpha,
	      BFmat,ot-1,umat,ot,&beta,BFtmpmat,ot-1);
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasTrans,ot-1,ot-1,ot,&alpha,
	      BFtmpmat,ot-1,BFmat,ot-1,&beta,retmat1,ot-1);
  delete [] BFtmpmat;
  delete [] BFmat;

  cout << "retmat0\n";
  printmat(retmat0,1,ot-1);

  cout << "retmat1\n";
  printmat(retmat1,ot-1,ot-1);

  cmplx** retarray=new cmplx*[2];
  retarray[0]=retmat0;
  retarray[1]=retmat1;
  return retarray;
}

cmplx* Xpsi(int ix, int ox, int ot, cmplx* Xmat, cmplx* psi, cmplx coeff){
  //psi is stored as a two dimensional array, with x as the first dimension
  //and t as the second.  This routine returns X.psi, where X is a spatial
  //matrix.
  //the relevant dimensions are
  //Xmat(ix,ox)
  //psi(ox,ot)
  //psiret(ix,ot)
  
  cmplx alpha=coeff;
  cmplx beta=0.;
  cmplx* psiret=new cmplx[ix*ot];
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ix,ot,ox,
	      &alpha,Xmat,ix,psi,ox,&beta,psiret,ix);
  return psiret;
  
}

cmplx* Tpsi(int it, int ot, int ox, cmplx* Tmat, cmplx* psi,cmplx coeff){
  //psi is stored as a two dimensional array, with x as the first dimension
  //and t as the second.  This routine returns psi.T, where T is a temporal
  //matrix  
  //the relevant dimensions are
  //Tmat(ot,it)
  //psi(ox,ot) 
  //psiret(ox,it)
  //cout << "ox "<<ox<<"\n";
  //cout << "ot "<<ot<<"\n";
  //cout << "it "<<it<<"\n";
  //cout << "psi\n";
  //printmat(psi,ox,ot);
  //cout << "Tpsi Tmat\n";
  //printmat(Tmat,ot,it);
  cmplx alpha=coeff;
  cmplx beta=0.;
  cmplx* psiret=new cmplx[ox*it];
  //cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ox,ot,it,
  //	      &alpha,psi,ox,Tmat,ot,&beta,psiret,ox);
  
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ox,it,ot,&alpha,
	      psi,ox,Tmat,ot,&beta,psiret,ox);
//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ox,it,ot,&alpha,psi,ox,Tmat,ot,&beta,psiret,ox);

  //cout << "Tpsi psiret\n";
  //printmat(psiret,ox,it);
  return psiret;
}

cmplx* XTpsi(int ix, int ox, int it, int ot, cmplx* Xmat, cmplx* Tmat, 
	     cmplx* psi,cmplx coeff){
  //multiply psi by Xmat and Tmat
  //relevant dimensions are
  //Xmat(ix,ox)
  //psi(ox,ot)
  //tmppsi(ix,ot)
  //Tmat(ot,it)
  //retpsi(ox,it)=Xmat.psi.Tmat
 
  //cout << "XTpsi ix,ox,it,ot "<<ix<<", "<<ox<<", "<<it<<", "<<ot<<"\n";

  cmplx* tmppsi=Xpsi(ix,ox,ot,Xmat,psi,coeff);
  //cout << "XTpsi tmppsi\n";
  //printmat(tmppsi,ix,ot);
  cmplx one=1.;
  //cout << "XTpsi Tmat\n";
  //printmat(Tmat,ot,it);

  cmplx* retpsi=Tpsi(it,ot,ox,Tmat,tmppsi,one);
  //cout << "XTpsi retpsi\n";
  //printmat(retpsi,ox,ot);
  delete [] tmppsi;
  return retpsi;
}

cmplx* Lpsi(int ox, int ot, cmplx** xmats, cmplx** tmats, 
	    cmplx* psi){
  //calculate (i*O.Q-T.U-A.p).psi = L.psi
  //(action residual = -L.psi)

  cout << "Lpsi input psi\n";
  printmat(psi,ox,ot);

  

  cmplx* umat=tmats[0];
  cmplx* qmat=tmats[1];
  cmplx* Amat=tmats[2];

  //cout << "Lpsi umat\n";
  //printmat(umat,ot,ot);
  //cout << "Lpsi qmat\n";
  //printmat(qmat,ot,ot);
  //cout << "Lpsi Amat\n";
  //printmat(Amat,ot,ot);

  cmplx* omat=xmats[0];
  cmplx* nablamat=xmats[1];
  cmplx* nablasqmat=xmats[2];

  

  cmplx OQcoeff=(0.,1.);
  //cmplx* __restrict OQpsi=XTpsi(ox,ox,ot,1,omat,qmat,psi,OQcoeff);
  cmplx* __restrict OQpsi=XTpsi(ox,ox,ot,ot,omat,qmat,psi,OQcoeff);

  cout << "Lpsi OQpsi\n";
  printmat(OQpsi,ox,ot);

  cmplx TUcoeff=-0.5;
  //cmplx* __restrict TUpsi=XTpsi(ox,ox,ot,1,nablasqmat,umat,psi,TUcoeff);
  cmplx* __restrict TUpsi=XTpsi(ox,ox,ot,ot,nablasqmat,umat,psi,TUcoeff);
  
  cout << "Lpsi TUpsi\n";
  printmat(TUpsi,ox,ot);

  cmplx APcoeff=-1.;//i A.p, where p=i*nabla
  //cmplx* __restrict APpsi=XTpsi(ox,ox,ot,1,nablamat,Amat,psi,APcoeff);
  cmplx* __restrict APpsi=XTpsi(ox,ox,ot,ot,nablamat,Amat,psi,APcoeff);

  cout << "Lpsi APpsi\n";
  printmat(APpsi,ox,ot);

  cmplx* __restrict psiret=initializearray(ox*ot);
  for(int i=0;i<ox*ot;i++){
    psiret[i]=OQpsi[i]+TUpsi[i]+APpsi[i];
  }

  //cout << "Lpsi psiret\n";
  //printmat(psiret,ox,ot);

  delete [] OQpsi;
  delete [] TUpsi;
  delete [] APpsi;
  return psiret;
}

cmplx* rhsvecsetup(int ot, int ox,cmplx** tmats, cmplx** xmats, cmplx* psi0){
  //set up the right hand side of the least action linear system (ie, the
  //residuals due to the initial condition psi0, projected onto the set of
  //temporal basis functions which are zero at the initial time).  
  //psi0 has shape (ox,1) 
  //spatial matrices have shape (ox,ox) 
  //temporal matrices have shape (1,(ot-1))
  cmplx* umatp=tmats[0];
  cmplx* qmatp=tmats[1];
  cmplx* Amatp=tmats[2];

  cout << "umatp\n";
  printmat(umatp,1,ot-1);
  cout << "qmatp\n";
  printmat(qmatp,1,ot-1);
  cout << "Amatp\n";
  printmat(Amatp,1,ot-1);

  cmplx* omat=xmats[0];
  cmplx* nablamat=xmats[1];
  cmplx* nablasqmat=xmats[2];

  cmplx OQcoeff=(0.,1.);
  cmplx* __restrict OQpsi=XTpsi(ox,ox,ot-1,1,omat,qmatp,psi0,OQcoeff);

  cmplx TUcoeff=-0.5;
  cmplx* __restrict TUpsi=XTpsi(ox,ox,ot-1,1,nablasqmat,umatp,psi0,TUcoeff);


  cmplx APcoeff=-1.;//i A.p, where p=i*nabla
  cmplx* __restrict APpsi=XTpsi(ox,ox,ot-1,1,nablamat,Amatp,psi0,APcoeff);


  //construct rhs vec, recalling that the linear system requires the t index
  //to correspond to the least significant index, while arrayindx makes t the
  //most significant index.

  cmplx* __restrict rhsvec=new cmplx[ox*(ot-1)];
  for(int i=0;i<ox;i++){
    for(int n=0;n<ot-1;n++){
      rhsvec[lsindx(i,n,ox,ot-1)]=OQpsi[arrayindx(i,n,ox,ot-1)]+
	TUpsi[arrayindx(i,n,ox,ot-1)]+APpsi[arrayindx(i,n,ox,ot-1)];
    }
  }
  delete [] OQpsi;
  delete [] TUpsi;
  delete [] APpsi;

  return rhsvec;
}

cmplx* XTmat(int ox1, int ox2, int ot1, int ot2, 
	     cmplx* __restrict Xmat,  cmplx* __restrict Tmat, cmplx coeff){
  //returns XTmat for use in setting up linear system, st.
  //XTmat[i,j,n,m]=Xmat[i,j]*Tmat[n,m]
  cmplx* __restrict retmat=new cmplx[ox1*ot1*ox2*ot2];
  
  for(int i=0;i<ox1;i++){
    for(int j=0;j<ox2;j++){
      cmplx xval=Xmat[arrayindx(i,j,ox1,ox2)];
      for(int n=0;n<ot1;n++){
	for(int m=0;m<ot2;m++){
	  cmplx tval=Tmat[arrayindx(n,m,ot1,ot2)];
	  int indx1=lsindx(i,n,ox1,ot1);
	  int indx2=lsindx(j,m,ox2,ot2);
	  int xtindx=arrayindx(indx1,indx2,ox1*ot1,ox2*ot2);
	  retmat[xtindx]=xval*tval*coeff;
	}
      }
    }
  }
  return retmat;
}

cmplx* lhsmatsetup(int otp, int ox,cmplx** tmats, cmplx** xmats){
  //set up the left hand side of the least action matrix.  Before calling this
  //toutine, tmats have been set up such that they have shape(otp,otp), where
  //otp=ot-1
  //spatial matrices have shape (ox,ox) 
  //temporal matrices have shape (otp,otp)
  cmplx* umatp=tmats[0];
  cmplx* qmatp=tmats[1];
  cmplx* Amatp=tmats[2];

  cmplx* omat=xmats[0];
  cmplx* nablamat=xmats[1];
  cmplx* nablasqmat=xmats[2];

  cmplx OQcoeff=(0.,1.);
  cmplx* __restrict OQmat=XTmat(ox,ox,otp,otp,omat,qmatp,OQcoeff);

  cmplx TUcoeff=-0.5;
  cmplx* __restrict TUmat=XTmat(ox,ox,otp,otp,nablasqmat,umatp,TUcoeff);

  cmplx APcoeff=-1.;//i A.p, where p=i*nabla
  cmplx* __restrict APmat=XTmat(ox,ox,otp,otp,nablamat,Amatp,APcoeff);

  int nbf=ox*otp;
  cmplx* __restrict lhsmat=new cmplx[nbf*nbf];
  for(int i=0;i<nbf*nbf;i++){
    lhsmat[i]=OQmat[i]+TUmat[i]+APmat[i];
  }
  delete [] OQmat;
  delete [] TUmat;
  delete [] APmat;
  return lhsmat;
}



cmplx** linearsystemsetup(derivedbasis* xbas,temporalbasis* tbas,cmplx* psi0){
  //given a basis in space and time, and an initial condition psi0, find the
  //lhs and rhs of the least action linear system.  Here, the zeroth temporal
  //basis function is assumed to be the only function which is nonzero at time
  //t0, and to be orthogonal to the other basis functions.  This will allow
  //solution of the relevant linear system without use of lagrange
  //multipliers, but will require projecting both sides of the linear system
  //onto the (nt-1) dimensional subspace of functions which are zero at t0
  //(ie, those coefficients which are not specified by the initial condition).

  //first, need to find temporal matrices which will be used in constructing
  //lhsmat, rhsvec
  
  int ot=tbas->order;
  int ox=xbas->order;

  cmplx* ICbfmat=leftborderfunctionbasis(ot);
  
  cmplx** ICumats=tmatsetup_ics(ot,tbas->umat,tbas->umat,ICbfmat);
  cmplx** ICqmats=tmatsetup_ics(ot,tbas->qmat,tbas->umat,ICbfmat);
  cmplx** ICAmats=tmatsetup_ics(ot,tbas->Amat,tbas->umat,ICbfmat);


  //container to pass x matrices to setup routines
  cmplx** xmats=new cmplx*[3];
  xmats[0]=xbas->omat;
  xmats[1]=xbas->nablamat;
  xmats[2]=xbas->nablasqmat;

  //the zeroth elements of ICumats, ICqmats, ICAmats will be used to construct
  //rhs vector
  cmplx** rhsvec_tmats=new cmplx*[3];
  rhsvec_tmats[0]=ICumats[0];
  rhsvec_tmats[1]=ICqmats[0];
  rhsvec_tmats[2]=ICAmats[0];

  cout << "rhsvec_tmats[0]\n";
  printmat(ICumats[0],1,ot-1);

  //the 1st elements of ICumats, ICqmats, ICamats will be used to construct
  //lhs matrix.
  cmplx** lhsmat_tmats=new cmplx*[3];
  lhsmat_tmats[0]=ICumats[1];
  lhsmat_tmats[1]=ICqmats[1];
  lhsmat_tmats[2]=ICAmats[1];

  cmplx* rhsvec=rhsvecsetup(ot,ox,rhsvec_tmats,xmats,psi0);
  cmplx* lhsmat=lhsmatsetup(ot-1,ox,lhsmat_tmats,xmats);
  
  cmplx** retarray=new cmplx*[2];
  retarray[0]=rhsvec;
  retarray[1]=lhsmat;

  delete [] xmats;
  delete [] rhsvec_tmats;
  delete [] lhsmat_tmats;
  delete [] ICbfmat;
  delete [] ICumats[0];
  delete [] ICumats[1];
  delete [] ICumats;
  delete [] ICqmats[0];
  delete [] ICqmats[1];
  delete [] ICqmats;
  delete [] ICAmats[0];
  delete [] ICAmats[1];
  delete [] ICAmats;
  

  return retarray;
}

