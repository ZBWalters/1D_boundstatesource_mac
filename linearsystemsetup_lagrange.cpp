#include "./classes.h"
#include "./linearsystemsetup.h"

//subroutines for setting up linear systems of equations using lagrange
//multipliers to specify initial conditions.

//The convention used here is that the index corresponding to the N+1st
//temporal basis fucntion corresponds to the lagrange multiplier for that
//spatial function.  Since we will be solving the dual problem, the initial
//condition for the wavefunction correction will by default be 0.

cmplx* lhsmatsetup_lagrange(int ot, int ox, cmplx** tmats, cmplx** xmats){
  //first set up L=i*O.Q-.5*T.U-A.p

  //Then use L to set up lhsmat for the linear system, which will differ from
  //L by including lagrange multipliers for every spatial basis function (so
  //that initial conditions may be specified).
  
  //read matrices from input arrays tmats and xmats  
  cmplx* umat=tmats[0];
  cmplx* qmat=tmats[1];
  cmplx* Amat=tmats[2];

  cout << "lhsmatsetup_lagrange tmats #2\n";
    cout << "umat\n";
    printmat(umat,ot,ot);
    cout << "qmat\n";
    printmat(qmat,ot,ot);
    cout << "Amat\n";
    printmat(Amat,ot,ot);

  cmplx* omat=xmats[0];
  cmplx* nablamat=xmats[1];
  cmplx* nablasqmat=xmats[2];
  
  //set up matrices for use in constructing least action linear system
  cmplx OQcoeff=(0.,1.);
  cmplx* __restrict OQmat=XTmat(ox,ox,ot,ot,omat,qmat,OQcoeff);

  cmplx TUcoeff=-0.5;
  cmplx* __restrict TUmat=XTmat(ox,ox,ot,ot,nablasqmat,umat,TUcoeff);

  cmplx APcoeff=-1.;//i A.p, where p=i*nabla
  cmplx* __restrict APmat=XTmat(ox,ox,ot,ot,nablamat,Amat,APcoeff);

  //set up linear system
  int neq=ox*(ot+1);
  cmplx* lhsmat=initializearray(neq*neq);


  //copy XT matrices into lhsmat
  for(int i=0;i<ox;i++){
    for(int n=0;n<ot;n++){
      int indx1=lsindx(i,n,ox,ot);
      int indx1p=lsindx(i,n,ox,ot+1);
      for(int j=0;j<ox;j++){
	for(int m=0;m<ot;m++){
	  int indx2=lsindx(j,m,ox,ot);
	  int indx2p=lsindx(j,m,ox,ot+1);
	  lhsmat[arrayindx(indx1p,indx2p,ox*(ot+1),ox*(ot+1))]=
	    OQmat[arrayindx(indx1,indx2,ox*ot,ox*ot)]+
	    TUmat[arrayindx(indx1,indx2,ox*ot,ox*ot)]+
	    APmat[arrayindx(indx1,indx2,ox*ot,ox*ot)];
	}
      }
    }
  }

 

  //set up lagrange multipliers
  for(int i=0;i<ox;i++){
    int m=ot;
    int indx2p=lsindx(i,m,ox,ot+1);
    for(int n=0;n<ot;n++){
      int indx1p=lsindx(i,n,ox,ot+1);
      lhsmat[arrayindx(indx1p,indx2p,ox*(ot+1),ox*(ot+1))]+=pow(-1.,n);
      lhsmat[arrayindx(indx2p,indx1p,ox*(ot+1),ox*(ot+1))]+=pow(-1.,n);
    }
  }

  cout << "lhsmat\n";
  //printmat_scipy(lhsmat,ox*(ot+1),ox*(ot+1));
  cout << "ox "<<ox<<"\n";
  cout << "ot "<<ot<<"\n";
  printmat(lhsmat,ox*(ot+1),ox*(ot+1));


  delete [] OQmat;
  delete [] TUmat;
  delete [] APmat;
  return lhsmat;
  
}

  
cmplx* rhsvecsetup_lagrange(int ot, int ox,cmplx* resid, cmplx* dpsi0){
  //solve for L.dpsi=-resid, with initial condition dpsi(t0)=dpsi0.  In most
  //cases, expect dpsi0=0, but since it doesn't cost anything extra, we might
  //as well program up the general case.

  cout << "rhsvecsetup_lagrange resid\n";
  printmat(resid,ox,ot);

  cmplx* rhsvec=initializearray(ox*(ot+1));
  for(int i=0;i<ox;i++){
    for(int n=0;n<ot;n++){
      int indx1=arrayindx(i,n,ox,ot);
      int indx1p=lsindx(i,n,ox,ot+1);
      rhsvec[indx1p]=resid[indx1];
    }
  }
  
  for(int i=0;i<ox;i++){
    int n=ot;
    int indx1p=lsindx(i,n,ox,ot+1);
    rhsvec[indx1p]=dpsi0[i];
  }

  return rhsvec;

}

cmplx* rhsvecsetup_lagrange(int ot, int ox, cmplx* resid){
  //same as rhsvecsetup_lagrange, but with implicit zero initial condition
  cmplx* psi0=initializearray(ox);
  cmplx* psiret=rhsvecsetup_lagrange(ot,ox,resid,psi0);
  delete [] psi0;
  return psiret;
}

cmplx** linearsystemsetup_lagrange(derivedbasis* xbas,temporalbasis* tbas,
				   cmplx* psi0){
  //set initial guess to be psi0*F0(t), with initial condition dpsi0=0 for the
  //correction.
  int ox=xbas->order;
  int ot=tbas->order;
  cmplx* psixt_guess=psi0xF0(ox,ot,ot,psi0);
  cmplx* dpsi0=initializearray(ox);
  
  cmplx** retarray=linearsystemsetup_lagrange(xbas,tbas,psixt_guess,dpsi0);
  delete [] psixt_guess;
  delete [] dpsi0;
  return retarray;
}

cmplx** linearsystemsetup_lagrange(derivedbasis* xbas,temporalbasis* tbas,
				   cmplx* psixt_guess,cmplx* dpsi0){
  //given a basis in space and time, and an initial condition psi0, find the
  //lhs and rhs of the least action linear system.  Here, the initial
  //condition is imposed using lagrange multipliers, so the temporal basis is
  //assumed to be a constant times a lagrange basis (so that \sum_{n} A_{i,n}
  //(-1)**n =psi0_{i}, where psi0=\sum_{i} psi0_{i} X_{i}(x)).

  //psi(x,t) may be found by first calculating the residual of psi=\sum_{i}
  //psi_{i} * f_0(t), where f_0(t) is the temporal basis funtion which has
  //nonzero value at time t0 & is orthogonal to all functions which are zero
  //at t0 (see notes "Estimating truncation error in the least action linear
  //system").  A correction to this initial guess is then found by solving for
  //L.dpsi=-r, with initial condition dpsi(t0)=0.  The advantage to doing
  //things this way is that you can also calculate the correction which would
  //result if f_0(t) were calculated to a lower order using the same method.
  //Thus, this approach allows you to estimate the effects of the last
  //temporal basis function, which may prove invaluable in calculating the
  //correct stepsize.

  int xorder=xbas->order;
  int torder=tbas->order;

  //container to pass x matrices to setup routines
  cmplx** xmats=new cmplx*[3];
  xmats[0]=xbas->omat;
  xmats[1]=xbas->nablamat;
  xmats[2]=xbas->nablasqmat;

  cmplx** tmats=new cmplx*[3];
  tmats[0]=tbas->umat;
  tmats[1]=tbas->qmat;
  tmats[2]=tbas->Amat;

  cout << "linearsystemsetup_lagrange xmats\n";
  printmat(xmats[0],xorder,xorder);
  printmat(xmats[1],xorder,xorder);
  printmat(xmats[2],xorder,xorder);

  cout << "linearsystemsetup_lagrange tmats\n";
  printmat(tmats[0],torder,torder);
  printmat(tmats[1],torder,torder);
  printmat(tmats[2],torder,torder);

  cout << "psixt_guess\n";
  printmat(psixt_guess,xorder,torder);

  //action residual of initial guess
  cmplx* resid=Lpsi(xorder,torder,xmats,tmats,psixt_guess);
  arrayscale(resid,xorder*torder,-1.);
  cout << "linearsystemsetup_lagrange resid\n";
  printmat(resid,xorder,torder);
    
  

//  //action residual due to initial condition
//  cmplx* psixt_ic=psi0xF0(xorder,torder,torder,psi0);
//  //action residual = -L.psi
//  cmplx* resid=Lpsi(xorder,torder,xmats,tmats,psixt_ic);
//  arrayscale(resid,xorder*torder,-1.);
//
//  cmplx* dpsi0=initializearray(xorder*torder);
  
  cmplx* rhsvec=rhsvecsetup_lagrange(torder,xorder,resid,dpsi0);
  cmplx* lhsmat=lhsmatsetup_lagrange(torder,xorder,tmats,xmats);

  cmplx** retarray=new cmplx*[2];
  retarray[0]=rhsvec;
  retarray[1]=lhsmat;

  //delete [] psixt_ic;
  //delete [] dpsi0;
  delete [] resid;
  //delete container arrays
  delete [] xmats;
  delete [] tmats;
  return retarray;

  

}

cmplx* psi0xF0(int xorder, int torder, int torderp,  cmplx* psi0){
  //when solving for a least action timestep, it is necessary to specify the
  //initial conditions.  For a certain order of temporal basis functions, this
  //can be done by constructing a function F0(t) which has value 1 at time t0
  //and is orthogonal to all polynomials of that order which are zero at time
  //t0.
  //This subroutine returns psi0(x)*F0(t), for a chosen order torderp of F0(t).
  

  //cout << "psi0xF0 xorder, torder, torderp "<<xorder<<", "<<torder<<", "<<torderp<<"\n";
  cmplx* basis=leftborderfunctionbasis(torderp);
  
  //cout << "left border function basis\n";
  //printmat(basis,torderp,torderp);

  cmplx* psiret=initializearray(xorder*torderp);

  //make psiret=psi0*F0(t)
  for(int i=0;i<xorder;i++){
    for(int n=0;n<torderp;n++){
      psiret[arrayindx(i,n,xorder,torder)]=
	psi0[i]*basis[arrayindx(0,n,torderp,torderp)];
    }
  }
  
  //cout << "psi0xF0 psiret\n";
  //printmat(psiret,xorder,torderp);

  delete [] basis;
  return psiret;
}
