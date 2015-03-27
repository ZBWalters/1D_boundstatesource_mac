#include "./classes.h"
#include "./wf.h"

void wf::propagate_iterative(pulse* Pls, potential* Pot,int steporder,rl tmin,rl tmax, rl dtfirst, rl accuracygoal){
  //propagate from tmin to tmax taking quality controlled steps
  rl t0=tmin;
  rl t1=t0;
  rl dt=dtfirst;
  rl newdt=dt;
  while(t0<tmax){
    t1=min(t0+dt,tmax);
    cout << "t0,t1,dt\t"<<t0<<"\t"<<t1<<"\t"<<dt<<"\n";
    gbas->updateHmats(Pot);
    timestep(t0,t1,steporder,accuracygoal,newdt);
    cout << "newdt after return\t"<<newdt<<"\n";
    t0=t1;
    dt=newdt;//min(newdt,.01);
  }
}

void wf::timestep(rl tmin, rl tmax,int steporder,rl accuracygoal,rl& newdt){
  //calculate wfxt which satisfies least action equations & calculate size of
  //next step
  rl dt=tmax-tmin;
  //update t mats
  //update Hamiltonians
  int nmax=50;

  int bfordergoal=steporder;
  int dtordergoal=steporder;//2.*steporder;
  cmplx* __restrict * __restrict * __restrict  steparrays=triD_iterativesolution(psi,bfordergoal,dtordergoal,nmax,dt);
  //printiterativesolution(nmax, steparrays[0],steparrays[1]);

  cmplx stepaction=iterativeaction(steparrays[0],steparrays[1],nmax);
  cout << "action "<<dt<<"\t"<<stepaction<<"\t"<<abs(stepaction)<<"\n";
  rl truncationerror=lastcomponentsignificance(steparrays[0],nmax);
  cout <<"last component significance\t"<<dt<<"\t"<<
    truncationerror<<"\n";

  //try new time step=dt*(accuracygoal/truncationerror)**(1./(bfordergoal-1))
  newdt=dt*pow(abs(accuracygoal/truncationerror),(1./(bfordergoal-1)));
  cout << "newdt\t"<<newdt<<"\n";
  delete [] psi;
  psi=psif(steparrays[0],nmax);
  
  //delete nonzero elements of guessxt, residxt
  for(int tmporder=0;tmporder<nmax;tmporder++){
    for(int i=0;i<2;i++){
      if(steparrays[i][tmporder]!=0){
	delete [] steparrays[i][tmporder];
      }
    }
  }
  delete [] steparrays[0];
  delete [] steparrays[1];
  delete [] steparrays;
}

cmplx* __restrict * __restrict * __restrict  wf::triD_iterativesolution(cmplx* __restrict  psi0, int bfordergoal, int dtordergoal,
				   int nmax, rl dt){
  //initialize guessxt, residxt
  //gbas->updateHmats();

  cmplx* __restrict * __restrict guessxt=new cmplx* __restrict [nmax];
  cmplx* __restrict * __restrict residxt=new cmplx* __restrict [nmax];
  for(int i=0;i<nmax;i++){
    guessxt[i]=0;//initializearray(nbf);
    residxt[i]=0;//initializearray(nbf);
  }

  //initial guess is that psixt=psi0 x F0
  triD_applydeltapsi(0,psi0,guessxt,residxt,dt);

  //for(int n=1;n<3;n++){
  //  correctioncycle(n,guessxt,residxt,dt);
  //}
  

  //initial guess has fractional error for all basis set coefficients.
  //correctbfsequence improves fractional error by a factor of dt for terms
  //from nmin to nmax.

//  //1st method:
//  //starting from n=0, work upwards to improve the guess until the fractional
//  //accuracy reaches dtorder (dtorder=1 from initial setup);
//  for (int dtorder=0;dtorder<dtordergoal;dtorder++){
//    correctbfsequence(1,bfordergoal,guessxt,residxt,dt);
//  }
//  //printiterativesolution(nmax,guessxt,residxt);

  //2nd method: 
  //starting from last nonzero residual, work downwards, so that
  //all residuals are of same order
  //cout << "dtordergoal "<<dtordergoal<<"\n";
  for(int lastindx=2;lastindx<dtordergoal;lastindx++){
    correctioncycle(lastindx,guessxt,residxt,dt);
    triD_updateresidual(lastindx,guessxt,residxt,dt);
  }

  cmplx* __restrict * __restrict * __restrict  retarray=new cmplx* __restrict * __restrict [2];
  retarray[0]=guessxt;
  retarray[1]=residxt;
  return retarray;
}

void wf::correctbfsequence(int nmin, int nmax, cmplx* __restrict * __restrict guessxt, 
			   cmplx* __restrict * __restrict residxt, rl dt){
  for(int n=nmin;n<=nmax;n++){
    cmplx* deltapsin=residualcorrection(n,guessxt,residxt,dt);
    //cout << "deltapsin from residualcorrection\n";
    //printmat(deltapsin,nbf,1);
    triD_applydeltapsi(n,deltapsin,guessxt,residxt,dt);
    delete [] deltapsin;
  }
}

void wf::correctioncycle(int n, cmplx* __restrict * __restrict guessxt, cmplx* __restrict * __restrict residxt, rl dt){
  //correct for the residual of order n, update guessxt and residxt, then do
  //the same thing for order n-1 until you've reached n=1
  
  //cout << "calling residualcorrection "<<n<<"\n";
  //cout <<"residxt[n]\n";
  //printmat(residxt[n],nbf,1);

  cmplx* deltapsin=residualcorrection(n,guessxt,residxt,dt);
  //cout << "deltapsin from residualcorrection\n";
  //printmat(deltapsin,nbf,1);
  triD_applydeltapsi(n,deltapsin,guessxt,residxt,dt);
  delete [] deltapsin;
  if(n>1){
    correctioncycle(n-1,guessxt,residxt,dt);
  }
}

cmplx* wf::residualcorrection(int n, cmplx* __restrict * __restrict guesst, cmplx* __restrict * __restrict residxt, rl dt){
  //set up and solve linear system for deltapsin s.t. 
  // (II Olap_{i,j} Q_{n,n}-H_{i,j} U_{n,n})C_{j,n}=-resid_{i,n}
  cmplx II=cmplx(0.,1.);
  int neq=nbf;

  //cout << "inside residualcorrection\n";

  //set up rhsmat to be -resid_{i,n}
  cmplx* __restrict rhsvec=new cmplx[neq];
  for(int i=0;i<neq;i++){
    rhsvec[i]=-residxt[n][i];
  }
  //cout << "rhsvec\n";
  //printmat(rhsvec,neq,1);
  //cout << "residxt[n]\n";
  //printmat(residxt[n],neq,1);

  int kl=gbas->xorder;
  int ku=gbas->xorder;
  int width=2*kl+ku+1;
  
  rl Qval=triD_Q(n,n,dt);
  rl Uval=triD_U(n,n,dt);

  //cout << "setting up lhsmat\n";
  //cout << "Qval, Uval "<<Qval<<", "<<Uval<<"\n";
  //cout << "neq "<<neq<<"\n";

  //copy element hmat, olap into lhsmat
  cmplx* __restrict lhsmat=new cmplx[nbf*width];
  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
    cmplx* __restrict tmph=tmpbasis->Hmat;
    cmplx* __restrict tmpo=tmpbasis->omat;
    int tmporder=tmpbasis->order;
    //cout << "tmph "<<eltnum<<"\n";
    //printmat(tmph,tmporder,tmporder);
    //cout << "tmpo "<<eltnum<<"\n";
    //printmat(tmpo,tmporder,tmporder);
    for(int i=0;i<tmporder;i++){
      for(int j=0;j<tmporder;j++){
	int indx1=gbas->eltfirstindices[eltnum]+i;
	int indx2=gbas->eltfirstindices[eltnum]+j;
	lhsmat[bandarrayindx(indx1,indx2,ku,kl,nbf)]+=
	  II*Qval*tmpo[arrayindx(i,j,tmporder,tmporder)]-
	  Uval*tmph[arrayindx(i,j,tmporder,tmporder)];
      }
    }
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

void wf::triD_updateresidual(int nmax,cmplx* __restrict * __restrict guessxt, 
			    cmplx* __restrict * __restrict residxt, rl dt){
  //set residxt to zero, guessxt to guessxt, residxt to L.guessxt.  this can
  //be done by first setting guessxt and residxt to zero, and residxt to
  //guessxt, then repeatedly calling triD_applydeltapsi to update guessxt and
  //residxt, before finally zeroing deltapsin (doing it this way allows re-use
  //of triD_applydeltapsi)
  
  cmplx * __restrict * __restrict deltapsin=new cmplx* __restrict[nmax];

  //first, move guessxt to deltapsin and set guessxt to zero
  for(int n=0;n<nmax;n++){
    //cout << "arraycopy "<<n<<"\n";
    deltapsin[n]=arraycopy(guessxt[n],nbf);
    zeroarray(guessxt[n],nbf);
  }
  //next, zero residxt (residxt for bf 0 is not needed)
  for(int n=1;n<=nmax;n++){
    //cout << "zero resid "<<n<<"\n";
    zeroarray(residxt[n],nbf);
  }
  //now, call triD_applydeltapsi to update guessxt and residxt
  for(int n=0;n<nmax;n++){
    //cout << "apply deltapsi"<<n<<"\n";
    triD_applydeltapsi(n,deltapsin[n],guessxt,residxt,dt);
    //finally, zero residxt[n];
    //cout << "delete deltapsi"<<n<<"\n";
    delete [] deltapsin[n];
  }
  delete [] deltapsin;
}


void wf::triD_applydeltapsi(int n, cmplx* __restrict deltapsin, 
			    cmplx* __restrict * __restrict guessxt, 
			    cmplx* __restrict * __restrict residxt, rl dt){
  //given a guessed wf guessxt, a correction deltapsin, and residuals residxt,
  //add deltapsin to guessxt, L.deltapsin to residxt
  cmplx II=cmplx(0.,1.);

  //apply update to guess
  if(guessxt[n]==0){
    guessxt[n]=initializearray(nbf);
  }
  for(int i=0;i<nbf;i++){
    guessxt[n][i]+=deltapsin[i];
  }
  


  //calculate residuals of deltapsin and apply to residxt
  cmplx* __restrict Hdeltapsi=gbas->Hpsi(deltapsin);
  cmplx* __restrict Olapdeltapsi=gbas->Olappsi(deltapsin);
  
  //cout << "Hdeltapsi\n";
  //printmat(Hdeltapsi,nbf,1);
  //
  //cout << "Olapdeltapsi\n";
  //printmat(Olapdeltapsi,nbf,1);

  
  int Uindxmin=max(1,n-1);
  int Uindxmax=n+1;
  int Qindxmin=1;
  int Qindxmax=max(1,n);
  //different definition of F0 means that n=0 is a special case
  if(n==0){
    Uindxmin=1;
    Uindxmax=2;
    Qindxmin=1;
    Qindxmax=2;
  }
  //apply II*Olap.Qmat
  for(int m=Qindxmin;m<=Qindxmax;m++){
    //cout << "n, m "<<n<<", "<<m<<", "<<Qindxmin<<", "<<Qindxmax<<"\n";
    rl Qval=triD_Q(n,m,dt);
    //cout << "Qval "<<Qval<<"\n";
    if(residxt[m]==0){
      residxt[m]=initializearray(nbf);
    }
    for(int j=0;j<nbf;j++){
      residxt[m][j]+=II*Qval*Olapdeltapsi[j];
      //cout << "delta j "<<II*Qval*Olapdeltapsi[j]<<"\n";
    }
    //cout << "residxt m "<<m<<"\n";
    //printmat(residxt[m],nbf,1);
  }
  //apply -H.Umat
  for(int m=Uindxmin;m<=Uindxmax;m++){
    rl Uval=triD_U(n,m,dt);
    if(residxt[m]==0){
      residxt[m]=initializearray(nbf);
    }
    for(int j=0;j<nbf;j++){
      residxt[m][j]-=Uval*Hdeltapsi[j];
    }
    //cout << "residxt m "<<m<<"\n";
    //printmat(residxt[m],nbf,1);
  }
  delete [] Hdeltapsi;
  delete [] Olapdeltapsi;
}

void wf::printiterativesolution(int nmax, cmplx* __restrict * __restrict guessxt, cmplx*__restrict * __restrict residxt){
  //print out current guess, residual, and ratio of residual to guess for all
  //orders where these are nonzero
  cout << "Printing current guess and residual\n";

  for(int tmporder=1;tmporder<nmax;tmporder++){
    if((guessxt[tmporder]!=0) or (residxt[tmporder]!=0)){
      cout << "order "<<tmporder<<"\n";
    }
    if(guessxt[tmporder]!=0){
      cout << "guessxt[order]\n";
      printmat(guessxt[tmporder],nbf,1);
    }
    if(residxt[tmporder]!=0){
      cout << "residxt[order]\n";
      printmat(residxt[tmporder],nbf,1);
    }
  }
}

cmplx wf::iterativeaction(cmplx* __restrict * __restrict guessxt,cmplx* __restrict * __restrict residxt, int nmax){
  //because residxt=L.guessxt, it's easy to calculate 
  //action =conj(guessxt).L.guessxt
  cmplx action=0.;
  for(int tmporder=0;tmporder<nmax;tmporder++){
    if((guessxt[tmporder]!=0) and (residxt[tmporder]!=0)){
      cmplx* __restrict tmpaction=new cmplx[1];
      cmplx alpha=1.;
      cmplx beta=0.;
      cblas_zgemv(CblasColMajor,CblasConjTrans,nbf,1,&alpha,guessxt[tmporder],
		  nbf,residxt[tmporder],1,&beta,tmpaction,1);
      action+=tmpaction[0];
      delete [] tmpaction;
    }
    
  }
  return action;
}

cmplx* wf::psif(cmplx* __restrict * __restrict guessxt, int nmax){
  //Since all but the zeroth of the tridiagonal temporal functions have value
  //1 at the end of the propagation step, psif is simply the sum of the
  //coefficients of guessxt
  //cout << "entered psif\t"<<nmax<<"\t"<<nbf<<"\n";
  cmplx* __restrict psiret=new cmplx[nbf];//initializearray(nbf);
  //cout << "initialized psiret\n";
  for(int tmporder=1;tmporder<nmax;tmporder++){
    //cout << "psif tmporder\t"<<tmporder<<"\n";
    if(guessxt[tmporder]!=0){
      for(int i=0;i<nbf;i++){
	psiret[i]+=guessxt[tmporder][i];
      }
    }
  }
  return psiret;
}

cmplx* wf::lastcomponentptr(cmplx* __restrict * __restrict xtarray, int nmax){
  //return last entry of guessxt, residxt which is not null
  int tmporder=nmax-1;
  bool loopflag=true;
  cmplx* retpointer=0;
  for(int i=0;i<nmax;i++){
    if(xtarray[i]!=0){
      retpointer=xtarray[i];
    }
  }
  return retpointer;
}

rl wf::lastcomponentsignificance(cmplx* __restrict * __restrict xtarray, int nmax){
  cmplx* __restrict wfsum=psif(xtarray,nmax);
  cmplx* __restrict lastcomponent=lastcomponentptr(xtarray,nmax);
  //rl retval= gbas->innerproduct(wfsum,lastcomponent);
  //for(int n=0;n<nmax;n++){
  //  cout << "significance of component "<<n<<"\t";
  //  cout <<sqrt(abs(gbas->norm(xtarray[n])))<<"\n";
  //}
  rl retval=sqrt(abs(gbas->norm(lastcomponent)));
  delete [] wfsum;
  return retval;
}

rl wf::nthcomponentsignificance(cmplx* __restrict * __restrict xtarray, int n){
  cmplx* __restrict nthcomponent=xtarray[n];
  //rl retval= gbas->innerproduct(wfsum,lastcomponent);
  //for(int n=0;n<nmax;n++){
  //  cout << "significance of component "<<n<<"\t";
  //  cout <<sqrt(abs(gbas->norm(xtarray[n])))<<"\n";
  //}
  rl retval=0.;
  if(nthcomponent!=0){
    retval=sqrt(abs(gbas->norm(nthcomponent)));
  }
  return retval;
}

rl wf::nthcomponentsignificance(cmplx* __restrict * __restrict xtarray1,
				cmplx* __restrict * __restrict xtarray2, 
				int n){
  cmplx* __restrict nthcomponent1=xtarray1[n];
  cmplx* __restrict nthcomponent2=xtarray2[n];
  cmplx* __restrict diff=initializearray(nbf);
  for(int i=0;i<nbf;i++){
    diff[i]=nthcomponent1[i]-nthcomponent2[i];
  }
  //rl retval= gbas->innerproduct(wfsum,lastcomponent);
  //for(int n=0;n<nmax;n++){
  //  cout << "significance of component "<<n<<"\t";
  //  cout <<sqrt(abs(gbas->norm(xtarray[n])))<<"\n";
  //}
  rl retval=sqrt(abs(gbas->norm(diff)));
  delete [] diff;
  return retval;
}

//rl wf::lastcomponentsignificance(cmplx* __restrict * __restrict xtarray, int nmax){
//  cmplx* __restrict wfsum=psif(xtarray,nmax);
//  cmplx* __restrict lastcomponent=lastcomponentptr(xtarray,nmax);
//  //rl retval= gbas->innerproduct(wfsum,lastcomponent);
//  rl retval=0.;
//  int nelts=gbas->nelts;
//  int* eltfirstindices=gbas->eltfirstindices;
//  for(int eltnum=0;eltnum<nelts;eltnum++){
//    cmplx* tmpsum=wfsum+eltfirstindices[eltnum];
//    cmplx* tmplast=lastcomponent+eltfirstindices[eltnum];
//    derivedbasis* tmpbas=gbas->elementbases[eltnum];
//    rl sumnorm=tmpbas->norm(tmpsum);
//    rl lastnorm=tmpbas->norm(tmplast);
//    //cout << "eltnum, lastnorm, sumnorm\t"<<eltnum<<"\t"<<lastnorm<<"\t"<<sumnorm<<"\n";
//    retval=max(retval,lastnorm/sumnorm);
//  }
//
//  //rl retval=sqrt(abs(gbas->norm(lastcomponent)));
//  delete [] wfsum;
//  return sqrt(retval);
//}



//rl wf::lastcomponentsignificance(cmplx* __restrict * __restrict xtarray, int nmax){
//  //find the spatial wavefunction corresponding to the last legendre
//  //polynomial, starting from a solution in the tridiagonal basis.  The
//  //difficulty with this approach is that the low order terms are so large
//  //relative to the high order terms that summing over all orders simply
//  //yields the errors in the low order terms, rather than anything useful.
//  int nlast=nmax-1;
//  cmplx coeff=0.;
//  int inc=1;
//  cmplx* __restrict lastcomponent=initializearray(nbf);
//  for(int n=0;n<nmax;n++){
//    if(n<2){
//      coeff=1.*pow(-1.,n+nlast);
//    }
//    else{
//      coeff=2.*pow(-1.,n+nlast);
//    }
//    cout <<"n,coeff "<<n<<", "<<coeff<<"\n";
//    //add coeff*xtarray[n] to lastcomponent
//    cblas_zaxpy(nbf,&coeff,xtarray[n],inc,lastcomponent,inc);
//    cout << "significance "<<(abs(gbas->norm(lastcomponent)))<<"\n";
//  }
//  rl retval=(abs(gbas->norm(lastcomponent)));
//  delete [] lastcomponent;
//  return retval;
//}
