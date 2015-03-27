#include "./classes.h"
#include "./wf.h"

//The routines in this file propagate the time dependent Schroedinger equation
//including an inhomogenous source term.  Such a source term arises when one
//wishes to calculate the deviation of the wavefunction from some
//almost-correct wavefunction -- in this context, the part of a wavefunction
//which tunnels free of a bound state.

void wf::propagate_iterative_inhomogenous(rl tmin, rl tmax, rl& dt,
					  pulse* Pls, potential* Pot,
					  wf* sourcewf, rl sourceEn, 
					  rl accuracygoal){
  //the iterative solution procedure works much easier in the tridiagonal
  //basis defined in triDbasis.cpp.  For each timestep set up this
  //temporalbasis.

  rl t0=tmin;
  rl t1=t0;
  rl newdt=dt;
  while(t0<tmax){
    t1=min(t0+dt,tmax);
    cout << "t0,t1,dt\t"<<t0<<"\t"<<t1<<"\t"<<dt<<"\n";
    //gbas->updateHmats(Pot);
    //sourcewf->gbas->updateHmats(Pot);
    //timestep_Uinv(t0,t1,Pot,Pls,sourcewf,sourceEn,accuracygoal,newdt);
    timestep_temporaleigenvectors(t0,t1,Pot,Pls,sourcewf,sourceEn,
				  accuracygoal,newdt);
    cout << "newdt after return\t"<<newdt<<"\n";
    t0=t1;
    dt=newdt;//min(newdt,.01);//test:maximum time step
  }
  
}

void wf::timestep(rl tmin, rl tmax, potential* Pot,pulse* Pls,
		  wf* sourcewf, rl sourceEn, rl accuracygoal, rl& newdt){
  //find wfxt which satisfies least action equations & calculate size of next
  //step

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
  //printmat(tbas_triD->umat,ot,ot);
  //cout << "qmat\n";
  //printmat(tbas_triD->qmat,ot,ot);
  //cout << "Amat\n";
  //printmat(tbas_triD->Amat,ot,ot);
  

  //targetxt is the action due to (delta H).sourcewf
  cmplx* __restrict * __restrict targetxt=sourcewf->
    inhomogeneous_action_lengthgauge(sourceEn,Pls,tbas_triD);
  //cmplx* __restrict * __restrict targetxt=sourcewf->
  //  inhomogeneous_action(sourceEn,tbas_triD);

  //cout << "targetxt\n";
  //for(int n=0;n<ot;n++){
  //  cout <<"n = "<<n<<"\n";
  //  printmat(targetxt[n],1,nbf);
  //}

  //initialize guessxt, residxt
  cmplx* __restrict * __restrict guessxt=new cmplx* __restrict[ot];
  cmplx* __restrict * __restrict residxt=new cmplx* __restrict[ot];
  for(int n=0;n<ot;n++){
    guessxt[n]=initializearray(nbf);
    residxt[n]=initializearray(nbf);
  }


  //initial guess is initial condition: guessxt[0]=psi
  applydeltapsi(0,ot,psi,guessxt,residxt,tbas_triD);
  //applydeltapsi(1,ot,psi,guessxt,residxt,tbas_triD);//choose initial guess
						    //that psi be constant




//  cout << "tbas_triD nablamat\n";
//  printmat(tbas_triD->nablamat,ot,ot);

//  //starting with n=1, apply correction cycle to iteratively improve the guess
//  correctionsequence(1,ot-1,ot,guessxt,residxt,targetxt,tbas_triD);
//  for(int n=1;n<ot;n++){
//    correctioncycle(n,ot,guessxt,residxt,targetxt,tbas_triD);
//  }
  
  //for(int iterindx=0;iterindx<2;iterindx++){

//  //cout << "n^2 operation\n";
//  for(int n=1;n<ot;n++){
//    //correctionsequence_dual(1,ot-1,ot,guessxt,residxt,targetxt,tbas_triD);
//    //correctionsequence_dual(ot-1,1,ot,guessxt,residxt,targetxt,tbas_triD);
//    correctionsequence_dual(n,1,ot,guessxt,residxt,targetxt,tbas_triD);
//  }

  

  cout << "refining\n";
  //initial condition of refinement cycle is that level 0 have 0 residual.
  //Then you adjust higher levels so that they have zero residuals without
  //paying attention to level 0.  End result will be a solution with zero
  //residual for nonzero levels, but no higher residual at level 0 than is
  //necessary to attain this
  correct_differentlevels(1,0,ot,guessxt,residxt,targetxt,tbas_triD);
  correctioncycle_dual_recursive(ot,1,ot,accuracygoal,guessxt,residxt,targetxt,
				 tbas_triD);
  //for(int iterindx=0;iterindx<5;iterindx++){
//  rl maxresid=1.;
//  int iterindx=0;
//  while(maxresid>accuracygoal){
//    iterindx++;
//    cout << "iterindx "<<iterindx<<"\n";
//    correctionsequence_dual(1,ot-1,ot,guessxt,residxt,targetxt,tbas_triD);
//    correctionsequence_dual(ot-1,1,ot,guessxt,residxt,targetxt,tbas_triD);
//    maxresid=0.;
//    for(int nt=1;nt<ot;nt++){
//      rl resid=nthcomponentsignificance(residxt,targetxt,nt);
//      maxresid=max(maxresid,resid);
//      cout << nt<<"\trefining residual significance\t"<<
//    	resid<<"\n";
//    }
//    cout << "maxresid\t"<<iterindx<<"\t"<<maxresid<<"\n";
//  }

//  int iterindx=0;
//  //while((errornorm(ot,residxt,targetxt) > accuracygoal) and (iterindx<ot)){
//  rl maxresidual=1.;
//  while((iterindx<ot) and (maxresidual>accuracygoal)){
//    iterindx++;
//    //  for(int iterindx=0;iterindx<2;iterindx++){
//    cout <<"iterindx "<<iterindx<<"\n";
//    int n=ot-1;
//    correctionsequence(1,ot-1,ot,guessxt,residxt,targetxt,tbas_triD);
//    correctioncycle(n,ot,guessxt,residxt,targetxt,tbas_triD);
//    //test agreement between residxt and targetxt
//    maxresidual=0.;
//    for(int nt=1;nt<ot;nt++){
//      rl significance=nthcomponentsignificance(residxt,targetxt,nt);
//      cout << nt<<" residual significance\t"<<
//	significance<<"\n";
//      maxresidual=max(maxresidual,significance);
//
//    }
//    cout << "max residual\t"<<maxresidual<<"\n";
//    rl enorm=errornorm(ot-1,residxt,targetxt);
//    cout <<"error norm\t"<<enorm<<"\n";
//    cout<<"\n";
//  }

  //cout << "guessxt\n";
  //for(int n=0;n<ot;n++){
  //  cout <<"n="<<n<<"\n";
  //  printmat(guessxt[n],nbf,1);
  //}

  rl enorm=errornorm(ot-1,residxt,targetxt);
  cout <<"error norm\t"<<enorm<<"\n";
  for(int nt=0;nt<ot;nt++){
    cout << nt<<"\tresidual significance\t"<<
      nthcomponentsignificance(residxt,targetxt,nt)<<"\n";
  }
  

  for(int nt=0;nt<ot;nt++){
    cout << nt<<"\tcomponent significance\t"<<
      nthcomponentsignificance(guessxt,nt)<<"\n";
  }

  

  rl truncationerror=lastcomponentsignificance(guessxt,ot);
  cout <<"last component significance\t"<<dt<<"\t"<<
    truncationerror<<"\n";

  //try new time step=dt*(accuracygoal/truncationerror)**(1./(bfordergoal-1))
  rl error=truncationerror;//max(truncationerror,maxresidual);
  newdt=dt*pow(abs(accuracygoal/error),(1./(ot-1)));
  cout << "newdt\t"<<newdt<<"\n";

  delete [] psi;
  psi=psif(guessxt,ot);

  
  delete [] tbas_triD_basischange;
  delete tbas_triD;
  for(int n=0;n<ot;n++){
    delete [] guessxt[n];
    delete [] residxt[n];
    delete [] targetxt[n];
  }
  delete [] guessxt;
  delete [] residxt;
  delete [] targetxt;
}

cmplx* wf::residualcorrection(int n, cmplx* __restrict resid, 
			      temporalbasis* tbas){
  return residualcorrection(n,n,resid,tbas);
}

cmplx* wf::residualcorrection(int ncorr,int nresid, cmplx* __restrict resid, 
			      temporalbasis* tbas){
  //solve for change to wf at order ncorr necessary to cancel a residual at
  //level nresid 
  //i,e, set up and solve linear system for deltapsin s.t.  
  //(II*Olap_{i,j} Q_{nresid,ncorr}-H_{i,j}
  //U_{nresid,ncorr}+A_{nresid,ncorr}.P_{i,ncorr}- 1/2
  //Asq_{nresid,ncorr}Olap_{i,j})C_{j,ncorr}=-resid_{i,n}
  
  cmplx II=cmplx(0.,1.);
  int neq=nbf;

  //cout << "inside residualcorrection\n";

  //set up rhsmat to be -resid_{i,n}
  cmplx* __restrict rhsvec=new cmplx[neq];
  for(int i=0;i<neq;i++){
    rhsvec[i]=-resid[i];
  }
  //cout << "rhsvec\n";
  //printmat(rhsvec,neq,1);
  //cout << "residxt[n]\n";
  //printmat(residxt[n],neq,1);


  int kl=gbas->xorder;
  int ku=gbas->xorder;
  int width=2*kl+ku+1;


  //find Amat.Umat
  int ot=tbas->order;
  cmplx* __restrict AUmat=initializearray(ot*ot);
  cmplx alpha=1.;
  cmplx beta=0.;
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
	      tbas->Amat,ot,tbas->umat,ot,&beta,AUmat,ot);
  //find Asqmat.Umat
  cmplx* __restrict AsqUmat=initializearray(ot*ot);
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
	      tbas->Amat,ot,AUmat,ot,&beta,AsqUmat,ot);
//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
//	      tbas->Asqmat,ot,tbas->umat,ot,&beta,AsqUmat,ot);

  cmplx Qval=tbas->qmat[arrayindx(nresid,ncorr,ot,ot)];
  cmplx Uval=tbas->umat[arrayindx(nresid,ncorr,ot,ot)];
  cmplx Aval=tbas->Amat[arrayindx(nresid,ncorr,ot,ot)];
  cmplx ATval=tbas->Amat[arrayindx(ncorr,nresid,ot,ot)];//value for A transpose
  cmplx Asqval=tbas->Asqmat[arrayindx(nresid,ncorr,ot,ot)];
  cmplx AUval=AUmat[arrayindx(nresid,ncorr,ot,ot)];
  cmplx AsqUval=AsqUmat[arrayindx(nresid,ncorr,ot,ot)];

  //if(ncorr!=nresid){
  //  cout << "ncorr\t"<<ncorr<<"\n";
  //  cout << "nresid\t"<<nresid<<"\n";
  //  cout << "Qval\t"<<Qval<<"\n";
  //  cout << "Uval\t"<<Uval<<"\n";
  //  cout << "Aval\t"<<Aval<<"\n";
  //  cout << "Asqval\t"<<Asqval<<"\n";
  //}

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
    
    //cout << "tmpnabla\n";
    //printmat(tmpnabla,tmporder,tmporder);
    //cout << "tmpnablasq\n";
    //printmat(tmpnablasq,tmporder,tmporder);

    //set up p.olap, olap.p for use in velocity gauge hamiltonian
    cmplx* __restrict tmpop=initializearray(tmporder*tmporder);
    cmplx* __restrict tmppo=initializearray(tmporder*tmporder);
    cmplx* __restrict tmposq=initializearray(tmporder*tmporder);

    cmplx alpha=cmplx(0.,1.);
    cmplx beta=0.;
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,tmporder,tmporder,
		tmporder,&alpha,tmpo,tmporder,tmpnabla,tmporder,&beta,
		tmpop,tmporder);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,tmporder,tmporder,
		tmporder,&alpha,tmpnabla,tmporder,tmpo,tmporder,&beta,
		tmppo,tmporder);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,tmporder,tmporder,
		tmporder,&alpha,tmpo,tmporder,tmpo,tmporder,&beta,
		tmposq,tmporder);

    
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
          Uval*(-.5*tmpnablasq[arrayindx(i,j,tmporder,tmporder)]+
		tmpV[arrayindx(i,j,tmporder,tmporder)])-
	  Aval*(II*tmpnabla[arrayindx(i,j,tmporder,tmporder)])+
	  Asqval*tmpo[arrayindx(i,j,tmporder,tmporder)];
        //lhsmat[bandarrayindx(indx1,indx2,ku,kl,nbf)]+=
        //  II*Qval*tmpo[arrayindx(i,j,tmporder,tmporder)]-
        //  Uval*(-.5*tmpnablasq[arrayindx(i,j,tmporder,tmporder)]+
	//	tmpV[arrayindx(i,j,tmporder,tmporder)])-
	//  0.5*ATval*II*tmpnabla[arrayindx(i,j,tmporder,tmporder)]+
	//  0.5*Aval*II*tmpnabla[arrayindx(j,i,tmporder,tmporder)]+
	//  Asqval*tmpo[arrayindx(i,j,tmporder,tmporder)];
	

	//AUval*0.5*(tmpop[arrayindx(i,j,tmporder,tmporder)]+
	  //	     tmppo[arrayindx(i,j,tmporder,tmporder)])-
	  //0.5*AsqUval*tmposq[arrayindx(i,j,tmporder,tmporder)];
      }
    }
    delete [] tmppo;
    delete [] tmpop;
    delete [] tmposq;
  }

  //solve for deltapsin using zgvsv_
  int info=0;
  int* ipiv=new int[neq];
  int nrhs=1;
  int ldab=2*kl+ku+1;

	zgbsv_(&neq, &kl, &ku, &nrhs, reinterpret_cast <Lcmplx*> (lhsmat), &ldab,
			ipiv, reinterpret_cast <Lcmplx*> (rhsvec), &neq, &info);
  if(info!=0){
    cout << "problem with zgvsv_ in residualcorrection! "<<info<<"\n";
  }
  delete [] AUmat;
  delete [] AsqUmat;
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

void wf::applydeltapsi(int nmax,cmplx* __restrict * __restrict deltapsi,
			 cmplx* __restrict * __restrict guessxt,
			 cmplx* __restrict * __restrict residxt, 
		       temporalbasis* tbas ){
  //call deltapsi for all temporal orders
  for(int n=0;n<nmax;n++){
    applydeltapsi(n,nmax,deltapsi[n],guessxt,residxt,tbas);
  }

}



void wf::correct_differentlevels(int ncorr, int nresid, int nmax,
				 cmplx* __restrict * __restrict  guessxt,
				 cmplx* __restrict * __restrict  residxt,
				 cmplx* __restrict * __restrict  targetxt,
				 temporalbasis* tbas){
  cmplx* __restrict diffxt=initializearray(nbf);
  for(int i=0;i<nbf;i++){
    diffxt[i]=(residxt[nresid][i]-targetxt[nresid][i]);
  }
  cmplx* deltapsi=residualcorrection(ncorr,nresid,diffxt,tbas);
  applydeltapsi(ncorr,nmax,deltapsi,guessxt,residxt,tbas);
  delete [] diffxt;
  delete [] deltapsi;
}

void wf::correction_leveln(int n,int nmax,cmplx* __restrict * __restrict  guessxt,
			 cmplx* __restrict * __restrict  residxt,
			 cmplx* __restrict * __restrict  targetxt,
			 temporalbasis* tbas){

  //find residual to be corrected for, find correction, apply it
  cmplx* __restrict diffxt=initializearray(nbf);
  for(int i=0;i<nbf;i++){
    diffxt[i]=(residxt[n][i]-targetxt[n][i]);
  }
  
  //cout << "n to be corrected "<<n<<"\n";
  //cout << "diffxt[n] before correction "<<gbas->norm(diffxt)<<"\n"; 
  
  cmplx* deltapsi=residualcorrection(n,diffxt,tbas);
  applydeltapsi(n,nmax,deltapsi,guessxt,residxt,tbas);

  //for(int i=0;i<nbf;i++){
  //  diffxt[i]=(residxt[n][i]-targetxt[n][i]);
  //}
  //cout << "diffxt[n] after correction "<<gbas->norm(diffxt)<<"\n";

  delete [] diffxt;
  delete [] deltapsi;
  
}

void wf::correctioncycle(int n,int nmax, cmplx* __restrict * __restrict  guessxt,
			 cmplx* __restrict * __restrict  residxt,
			 cmplx* __restrict * __restrict  targetxt,
			 temporalbasis* tbas){
  //given a target action residual targetxt, apply corrections to temporal
  //order n, then recursively call to correct order n-1, etc.  Since n=0
  //corresponds to the initial condition for the tridiagonal temporal basis,
  //this level is exact and should not be altered.  For every n, calculate the
  //correction assuming the temporal matrices are diagonal -- this will not
  //solve exactly, but will yield a rapidly decreasing error term

  //first approach: correct guess at level n to cancel residuals at level n
  //cmplx* __restrict diffxt=initializearray(nbf);
  //for(int i=0;i<nbf;i++){
  //  diffxt[i]=(residxt[n][i]-targetxt[n][i]);
  //}
  correction_leveln(n,nmax,guessxt,residxt,targetxt,tbas);

//  //second approach: correct guess at level n to cancel residuals at level n-1
//  cmplx* __restrict diffxt=initializearray(nbf);
//  for(int i=0;i<nbf;i++){
//    diffxt[i]=(residxt[n-1][i]-targetxt[n][i]);
//  }
//  correct_differentlevels(n,n-1,nmax,guessxt,residxt,targetxt,tbas);

  //delete [] diffxt;
  if(n>1){
    correctioncycle(n-1,nmax,guessxt,residxt,targetxt,tbas);
  }
}

void wf::correctioncycle_dual_recursive(int nmax,int recursionlevel, 
					int maxrecursionlevel,rl accuracygoal,
					cmplx* __restrict * __restrict guessxt,
					cmplx* __restrict * __restrict residxt,
					cmplx* __restrict * __restrict targetxt,
					temporalbasis* tbas){
  //set up difference between guess & target & apply correction cycle.  If solution doesn't meet accuracy goal, recursively call self
  int ot=tbas->order;
  cmplx* __restrict * __restrict deltaxt=new cmplx* __restrict[ot];
  cmplx* __restrict * __restrict dresid=new cmplx* __restrict[ot];
  cmplx* __restrict * __restrict target2=new cmplx* __restrict[ot];
  for(int n=0;n<ot;n++){
    deltaxt[n]=initializearray(nbf);
    dresid[n]=initializearray(nbf);
    target2[n]=initializearray(nbf);
    for(int i=0;i<nbf;i++){
      target2[n][i]=targetxt[n][i]-residxt[n][i];
    }
    //cout << n<<"\ttarget significance\t"<<
    //  nthcomponentsignificance(targetxt,n)<<"\n";
  }
  

  ////cancel residuals from level 1 to ot-1, ascending
  //correctionsequence_dual(1,recursionlevel,ot,deltaxt,dresid,target2,tbas);
  correctionsequence_dual(1,ot-1,ot,deltaxt,dresid,target2,tbas);
  ////cancel residuals from level ot-1 to 1, descending
  correctionsequence_dual(ot-2,1,ot,deltaxt,dresid,target2,tbas);
  //correctionsequence_dual(recursionlevel,1,ot,deltaxt,dresid,target2,tbas);
  

  rl maxdelta=0.;
  for(int n=1;n<ot;n++){
    //rl delta=nthcomponentsignificance(deltaxt,n);
    rl delta=nthcomponentsignificance(dresid,target2,n);
    maxdelta=max(maxdelta,delta);
    cout << n<<"\ttarget-dresid significance\t"<<recursionlevel<<"\t"<<
      delta<<"\n";
  }
  cout << "maxdelta\t"<<maxdelta<<"\n";
  for(int n=1;n<ot;n++){
    rl delta=nthcomponentsignificance(deltaxt,n);
    cout << n<<"\tdeltaxt significance\t"<<recursionlevel<<"\t"<<
      delta<<"\n";
  }
  

  //if accuracy not good enough, recursively call self
  if((maxdelta>accuracygoal) and ((recursionlevel+1)<maxrecursionlevel)){
    correctioncycle_dual_recursive(nmax,recursionlevel+1,maxrecursionlevel,
				   accuracygoal,deltaxt,dresid,
				   target2,tbas);
  }

  //add deltaxt to guessxt
  for(int n=0;n<ot;n++){
    for(int i=0;i<nbf;i++){
      guessxt[n][i]+=deltaxt[n][i];
      residxt[n][i]+=dresid[n][i];
    }
  }

  //delete temporary arrays
  for(int n=0;n<ot;n++){
    delete [] deltaxt[n];
    delete [] dresid[n];
    delete [] target2[n];
  }
  delete [] deltaxt;
  delete [] dresid;
  delete [] target2;
  
  
}


void wf::correctionsequence_dual(int n1, int n2,int nmax, 
				 cmplx* __restrict * __restrict  guessxt,
				 cmplx* __restrict * __restrict  residxt,
				 cmplx* __restrict * __restrict  targetxt,
				 temporalbasis* tbas){
  //set up dual problem:
  int ot=tbas->order;
  cmplx* __restrict * __restrict deltaxt=new cmplx* __restrict[ot];
  cmplx* __restrict * __restrict dresid=new cmplx* __restrict[ot];
  for(int n=0;n<ot;n++){
    deltaxt[n]=initializearray(nbf);
    dresid[n]=initializearray(nbf);
  }
  cmplx* __restrict * __restrict target2=residualdiff(guessxt,
						      targetxt,tbas);
  correctionsequence(n1,n2,nmax,deltaxt,dresid,target2,tbas);
  applydeltapsi(nmax,deltaxt,guessxt,residxt,tbas);
  
  //for(int n=0;n<ot;n++){
  //  cout << "deltaxt "<<n<<" norm\t"<<nthcomponentsignificance(deltaxt,n)
  //	 <<"\t"<<nthcomponentsignificance(guessxt,n)<<"\n";
  //}

  //delete deltaxt,dresid,target2
  for(int n=0;n<ot;n++){
    delete [] deltaxt[n];
    delete [] dresid[n];
    delete [] target2[n];
  }
}

void wf::correctionsequence(int n1, int n2,int nmax, 
			    cmplx* __restrict * __restrict  guessxt,
			    cmplx* __restrict * __restrict  residxt,
			    cmplx* __restrict * __restrict  targetxt,
			    temporalbasis* tbas){
  //apply corrections to levels n1 to n2, where n2>=n1.  This yields bottom-up
  //correction, as opposed to the top-down correction used in correctioncycle.
  //Bottom-up correction is useful for improving a solution by increasing the
  //order of every residual by 1.
  if(n1<=n2){
    for(int n=n1;n<=n2;n++){
      //correct_differentlevels(n,n-1,nmax,guessxt,residxt,targetxt,tbas);
      correct_differentlevels(n,n,nmax,guessxt,residxt,targetxt,tbas);
      //correction_leveln(n,nmax,guessxt,residxt,targetxt,tbas);
    }
  }
  if(n1>n2){
    for(int n=n1;n>=n2;n--){
      //correct_differentlevels(n,n-1,nmax,guessxt,residxt,targetxt,tbas);
      correct_differentlevels(n,n,nmax,guessxt,residxt,targetxt,tbas);
      //correction_leveln(n,nmax,guessxt,residxt,targetxt,tbas);
    }
  }
  
}

cmplx wf::iterativeaction_inhomogenous(int nmax,
				       cmplx* __restrict * __restrict  guessxt,
				       cmplx* __restrict * __restrict  residxt,
				       cmplx* __restrict * __restrict  targetxt)
{
  //The error due to propagation with an inhomogenous source term phi=(delta
  //H).psi0(t) is given by (delta psi)^(HC).(L.delta psi)-(delta psi)^(HC).phi
  //Here L.delta psi is residxt, while phi is targetxt, so the action integral
  //is action=(delta psi)^(HC).(residxt-targetxt)

  cmplx alpha=1.;
  cmplx beta=0.;
  cmplx* __restrict actionarray=new cmplx[nmax];//store actions for later
						//summation -- allows loop to
						//be parallelized.
  for(int n=0;n<nmax;n++){
    cmplx* __restrict tmpaction=new cmplx[1];
    cmplx* __restrict tmpdiff=new cmplx[nbf];
    for(int i=0;i<nbf;i++){
      tmpdiff[i]=residxt[n][i]-targetxt[n][i];
    }
    cblas_zgemv(CblasColMajor,CblasConjTrans,nbf,1,&alpha,guessxt[n],
                  nbf,tmpdiff,1,&beta,tmpaction,1);
    actionarray[n]=tmpaction[0];
    delete [] tmpdiff;
    delete [] tmpaction;
  }

  //sum over actionarray to get full action
  cmplx action=0.;
  for(int n=0;n<nmax;n++){
    action+=actionarray[n];
  }
  delete [] actionarray;

  return action;
}

rl wf::errornorm(int nmax,cmplx* __restrict * __restrict  residxt,
	       cmplx* __restrict * __restrict  targetxt){
  rl retval=0.;
  for(int n=1;n<nmax;n++){
    cmplx* tmpdiff=arraydiff(targetxt[n],residxt[n],nbf);
    rl tmpnorm=gbas->norm(tmpdiff);
    retval+=tmpnorm;
    //cout << "error norm "<<n<<"\t"<<tmpnorm<<"\n";
    //cout << "target norm "<<gbas->norm(targetxt[n])<<"\n";
    delete [] tmpdiff;
  }
  //cout <<"\n";
  return retval;
}

cmplx* __restrict * __restrict 
wf::residualdiff(cmplx* __restrict * __restrict guessxt,
		 cmplx* __restrict * __restrict targetxt,temporalbasis* tbas){

  int ot=tbas->order;

  //initialize guessxt, residxt
  cmplx* __restrict * __restrict wfxt=new cmplx* __restrict[ot];
  cmplx* __restrict * __restrict rhsxt=new cmplx* __restrict[ot];
  cmplx* __restrict * __restrict retarray=new cmplx* __restrict[ot];

  for(int n=0;n<ot;n++){
    wfxt[n]=initializearray(nbf);
    rhsxt[n]=initializearray(nbf);
    retarray[n]=initializearray(nbf);
  }


  //find residual corresponding to guessxt
  for(int n=0;n<ot;n++){
    applydeltapsi(n,ot,guessxt[n],wfxt,rhsxt,tbas);
  }

  //return difference between residual and target
  for(int n=0;n<ot;n++){
    for(int i=0;i<nbf;i++){
      retarray[n][i]=targetxt[n][i]-rhsxt[n][i];
    }
  }

  //delete wfxt,rhsxt
  for(int n=0;n<ot;n++){
    delete [] wfxt[n];
    delete [] rhsxt[n];
  }
  delete [] wfxt;
  delete [] rhsxt;
  
  return retarray;
}

