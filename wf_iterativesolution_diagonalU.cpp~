#include "./classes.h"
#include "./wf.h"


//This file contains routines designed to iteratively solve the least action
//equation (i O Q - H U) (dpsi)=R by inserting Uinv.U between the lagrangian
//and the wavefunction, so that the new equation becomes 
//(i O (Q Uinv) - H (U Uinv)) (U dpsi)=R

//The new equation has the property that U Uinv is diagonal in the time index,
//while Q Uinv is lower triangular.  Because of this, the correction needed to
//zero residuals for time indices greater than or equal to 1 may be found by
//backsubstitution, so that the overall iterative procedure requires only n
//linear solves, as opposed to n^2 for previous way.  In addition, careful
//choice of the initial guess for dpsi can make it so that the residual of the
//initial guess does not include any terms of the form H*psi, so that the
//overall iterative procedure will not involve repeated multiplications by H.
//This will help numerical stability for high energy eigenvalues.

void wf::Upsi0(int nmax, cmplx* __restrict psi0, 
		cmplx* __restrict * __restrict psixt,
		cmplx* __restrict * __restrict residxt,rl scale,
		temporalbasis* tbas){
	//cmplx* vec=ICvec_(nmax);
	int ot=tbas->order;
	rl dt=tbas->tmax-tbas->tmin;
	cmplx* vec=initializearray(nmax);
	vec[0]=1.;
	vec[1]=1;//choose initial guess to be constant function
	for(int n=0;n<nmax;n++){
		zeroarray(psixt[n],nbf);
		zeroarray(residxt[n],nbf);
	}

	//read initial guess into tmppsixt
	cmplx* __restrict * __restrict tmppsixt=new cmplx*[ot];
	tmppsixt[0]=arraycopy(psi0,nbf);
	tmppsixt[1]=arraycopy(psi0,nbf);
	for(int n=2;n<ot;n++){
		tmppsixt[n]=initializearray(nbf);
	}


	cmplx* __restrict * __restrict Utmppsixt=Upsi(ot,tmppsixt,scale);
	//now apply Upsi to guessxt and resid
	for(int n=0;n<ot;n++){
		applydeltapsi_Uinv(n,ot, Utmppsixt[n],psixt,residxt,scale,tbas);
	}

	//delete temporary arrays
	for(int n=0;n<ot;n++){
		delete [] tmppsixt[n];
		delete [] Utmppsixt[n];
	}
	delete [] tmppsixt;
	delete [] Utmppsixt;


}
void wf::initialguess_Uinv(int nmax, cmplx* __restrict psi0, 
		cmplx* __restrict * __restrict psixt,
		cmplx* __restrict * __restrict resid,
		temporalbasis* tbas){
	rl scale=1.;
	initialguess_Uinv(nmax,psi0,psixt,resid,scale,tbas);
}


void wf::initialguess_Uinv(int nmax, cmplx* __restrict psi0, 
		cmplx* __restrict * __restrict psixt,
		cmplx* __restrict * __restrict resid,rl scale,
		temporalbasis* tbas){
	//choose psixt s.t. psixt[0]=psi0, Sum_m Unm.psixt[m]=0 for n>0.  This will
	//help with numerical stability, since it will mean that H.U.psixt=0 for
	//n>0, so that there's no repeated multiplication by Hamiltonian matrix to
	//worry about.

	//cmplx* vec=ICvec_(nmax);
	rl dt=tbas->tmax-tbas->tmin;
	cmplx* vec=initializearray(nmax);
	vec[0]=1.;
	vec[1]=1;//choose initial guess to be constant function
	for(int n=0;n<nmax;n++){
		zeroarray(psixt[n],nbf);
		zeroarray(resid[n],nbf);
	}

	for(int n=0;n<nmax;n++){
		cmplx* __restrict tmppsi=new cmplx[nbf];
		for(int i=0;i<nbf;i++){
			tmppsi[i]=psi0[i]*vec[n];
		}
		applydeltapsi_Uinv(n,nmax,tmppsi,psixt,resid,scale,tbas);
	}

}

cmplx* __restrict * __restrict  
wf::Uinvpsi(int nmax,cmplx* __restrict * __restrict dphi){
	//Given dphi=U.dpsix, solve for dpsixt
	cmplx* __restrict * __restrict dpsixt=new cmplx* __restrict[nmax];
	for(int n=0;n<nmax;n++){
		dpsixt[n]=initializearray(nbf);
	}
	cmplx* __restrict rhsmat=new cmplx[nmax*nbf];
	for(int n=0;n<nmax;n++){
		for(int i=0;i<nbf;i++){
			//rhsmat[arrayindx(i,n,nbf,nmax)]=dphi[n][i];
			rhsmat[arrayindx(n,i,nmax,nbf)]=dphi[n][i];
		}
	}


	//for every spatial index, solve U_nm.dpsixt[m][i]=dphi[n][i]

	//take advantage of tridiagonal structure to solve using gaussian elimination
	cmplx* __restrict Umat=triD_Umat(nmax);
	cmplx* __restrict Diag=new cmplx[nmax];
	cmplx* __restrict DL=new cmplx[nmax-1];
	cmplx* __restrict DU=new cmplx[nmax-1];
	for(int n=0;n<nmax;n++){
		Diag[n]=triD_U(n,n);//Umat[arrayindx(n,n,nmax,nmax)];
	}
	for(int n=1;n<nmax;n++){
		DU[n-1]=triD_U(n,n-1);//Umat[arrayindx(n,n-1,nmax,nmax)];
		DL[n-1]=triD_U(n-1,n);//Umat[arrayindx(n-1,n,nmax,nmax)];
	}
	int nrhs=nbf;
	int info=0;
	zgtsv_(&nmax, &nrhs, reinterpret_cast <Lcmplx*> (DL), reinterpret_cast
			<Lcmplx*> (Diag), reinterpret_cast <Lcmplx*> (DU), reinterpret_cast
			<Lcmplx*> (rhsmat), &nmax, &info);
	if(info!=0){
		cout << "problem with zgtsv_!\n";
	}
	//cout << "back from zgtsv_\n";





	//  //LU factorization method
	//  //LU factorization of Umat will be used by all solvers
	//  cmplx* __restrict UmatLU=triD_Umat(nmax);
	//  int info=0;
	//  int* ipiv=new int[nmax];
	//  zgetrf(&nmax,&nmax,UmatLU,&nmax,ipiv,&info);
	//  if(info!=0){
	//    cout << "problem with zgetrf in Uinvpsi!\n";
	//  }
	//    //solve for dpsixt using backsubstitution
	//  int nrhs=nbf;
	//  char trans='N';
	//  zgetrs(&trans,&nmax,&nrhs,UmatLU,&nmax,ipiv,rhsmat,&nmax,&info);
	//  if(info!=0){
	//    cout << "problem with zgetrs in Uinvpsi!\n";
	//  }

	for(int n=0;n<nmax;n++){
		for(int i=0;i<nbf;i++){
			//dpsixt[n][i]=rhsmat[arrayindx(i,n,nbf,nmax)];
			dpsixt[n][i]=rhsmat[arrayindx(n,i,nmax,nbf)];
		}
	}

	//cout << "deleting\n";
	delete [] rhsmat;
	//cout << "rhsmat\n";
	//delete [] UmatLU;
	//delete [] ipiv;
	delete [] Diag;
	//cout << "Diag\n";
	delete [] DL;
	//cout << "DL\n";
	delete [] DU;
	delete [] Umat;


	return dpsixt;
}

cmplx* __restrict * __restrict  
wf::Qinvpsi(int nmax,cmplx* __restrict * __restrict dphi, rl scale){
	//find Q^(-1)*dphi in tridiagonal basis assuming that dphi has been scaled
	//by scalefactor, so that dphi[n] has a prefactor of scale**-n.  Note that Q
	//is not invertible if order zero is included (the derivative of a constant
	//is zero), but is invertible in the tridiagonal basis using the range of
	//orders [1,infinity).  In this restricted range, Q^-1 has a simple upper
	//triangular form, with 2 on the diagonal and 4*(-1)**(m-n) for m>n.  Thus,
	//this routine doesn't do a matrix inversion, but just applies this solution
	//to dphi


	//initialize and zero retmat
	cmplx* __restrict * __restrict retmat=new cmplx*[nmax];
	for(int n=0;n<nmax;n++){
		retmat[n]=initializearray(nbf);
	}

	//for n>0, (Qinv.psi)[n][i]=2 psi[n][i]+ 
	//Sum_{m>n} 4 (-1)**(m-n) scale**(m-n) psi[m][i]
	cmplx* __restrict summat=initializearray(nbf);
	for(int n=(nmax-1);n>0;n--){
		for(int i=0;i<nbf;i++){
			retmat[n][i]=2.*dphi[n][i]-4.*scale*summat[i];
		}
		for(int i=0;i<nbf;i++){
			summat[i]=dphi[n][i]-scale*summat[i];
		} 
	}
	delete [] summat;
	return retmat;

}

cmplx* __restrict * __restrict  
wf::Upsi(int nmax,cmplx* __restrict * __restrict dphi, rl scale){
	//find U.psi, where U does not include the factor of dt

	//initialize and zero retmat
	cmplx* __restrict * __restrict retmat=new cmplx*[nmax];
	for(int n=0;n<nmax;n++){
		retmat[n]=initializearray(nbf);
	}

	//sub diagonal
	for(int n=1;n<nmax;n++){
		for(int i=0;i<nbf;i++){
			retmat[n][i]+=dphi[n-1][i]*triD_U(n,n-1)/scale;
		}
	}
	//main diagonal
	for(int n=0;n<nmax;n++){
		for(int i=0;i<nbf;i++){
			retmat[n][i]+=dphi[n][i]*triD_U(n,n);
		}
	}
	//super diagonal
	for(int n=0;n<nmax-1;n++){
		//cout << "superdiag n\t"<<n<<"\t"<<triD_U(n,n+1)<<"\n";
		for(int i=0;i<nbf;i++){
			retmat[n][i]+=dphi[n+1][i]*triD_U(n,n+1)/scale;
		}
	}
	return retmat;
}

cmplx* wf::residualcorrection_Uinv(int ncorr,int nresid, cmplx* __restrict resid, 
		temporalbasis* tbas){
	//solve for change to wf at order ncorr necessary to cancel a residual at
	//level nresid 
	//i,e, set up and solve linear system for deltapsin s.t.  
	//((II*Olap_{i,j} Q_{nresid,ncorr}-H_{i,j}
	//U_{nresid,ncorr}+A_{nresid,ncorr}.P_{i,ncorr}- 1/2
	//Asq_{nresid,ncorr}Olap_{i,j}) Uinv_jk) (U_kl C_{l,ncorr})=-resid_{i,n}

	//cout << "inside residualcorrection_Uinv\n";

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

	//find Uinv for tridiagonal basis
	int ot=tbas->order;
	//cout << "calling triD_Uinv\n";
	cmplx* __restrict Uinv=triD_Uinv(ot);
	cmplx alpha=1.;
	cmplx beta=0.;
	//cout << "found Uinv\n";

	//find Uinv.Qmat
	cmplx* __restrict UinvQ=initializearray(ot*ot);
	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
			Uinv,ot,tbas->qmat,ot,&beta,UinvQ,ot);
	//cout << "found QUinv\n";

	//find Uinv.Umat(should be diagonal w/factor dt on diagonals)
	cmplx* __restrict UinvU=initializearray(ot*ot);
	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
			Uinv,ot,tbas->umat,ot,&beta,UinvU,ot);
	//cout << "found UUinv\n";

	//find Uinv.Amat  
	cmplx* __restrict UinvA=initializearray(ot*ot);
	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
			Uinv,ot,tbas->Amat,ot,&beta,UinvA,ot);
	//cout << "found AUinv\n";

	//find Uinv.Asqmat
	cmplx* __restrict UinvAsq=initializearray(ot*ot);
	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
			Uinv,ot,tbas->Asqmat,ot,&beta,UinvAsq,ot);
	//cout << "found AsqUinv\n";
	//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
	//	      tbas->Asqmat,ot,tbas->umat,ot,&beta,AsqUmat,ot);

	cmplx Qval=UinvQ[arrayindx(nresid,ncorr,ot,ot)];
	cmplx Uval=UinvU[arrayindx(nresid,ncorr,ot,ot)];
	cmplx Aval=UinvA[arrayindx(nresid,ncorr,ot,ot)];
	cmplx Asqval=UinvAsq[arrayindx(nresid,ncorr,ot,ot)];
	//cmplx ATval=tbas->Amat[arrayindx(ncorr,nresid,ot,ot)];//value for A transpose

	//cmplx AUval=AUmat[arrayindx(nresid,ncorr,ot,ot)];
	//cmplx AsqUval=AsqUmat[arrayindx(nresid,ncorr,ot,ot)];

	//if(ncorr!=nresid){
	cout << "ncorr\t"<<ncorr<<"\n";
	cout << "nresid\t"<<nresid<<"\n";
	cout << "Qval\t"<<Qval<<"\n";
	cout << "Uval\t"<<Uval<<"\n";
	cout << "Aval\t"<<Aval<<"\n";
	cout << "Asqval\t"<<Asqval<<"\n";
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
			}
		}
		delete [] tmppo;
		delete [] tmpop;
		delete [] tmposq;
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
	delete [] Uinv;
	delete [] UinvQ;
	delete [] UinvU;
	delete [] UinvA;
	delete [] UinvAsq;
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

void wf::applydeltapsi_Uinv(int nmax,cmplx* __restrict * __restrict deltapsi,
		cmplx* __restrict * __restrict guessxt,
		cmplx* __restrict * __restrict residxt,
		temporalbasis* tbas ){
	//call deltapsi for all temporal orders
	rl scale=1.;
	for(int n=0;n<nmax;n++){
		applydeltapsi_Uinv(n,nmax,deltapsi[n],guessxt,residxt,scale,tbas);
	}

}

void wf::applydeltapsi_Uinv(int nmax,cmplx* __restrict * __restrict deltapsi,
		cmplx* __restrict * __restrict guessxt,
		cmplx* __restrict * __restrict residxt, rl scale,
		temporalbasis* tbas ){
	//call deltapsi for all temporal orders
	for(int n=0;n<nmax;n++){
		applydeltapsi_Uinv(n,nmax,deltapsi[n],guessxt,residxt,scale,tbas);
	}

}

void wf::applydeltapsi_Uinv(int ncorr,int nmax, cmplx* __restrict deltapsi,
		cmplx* __restrict * __restrict guessxt,
		cmplx* __restrict * __restrict residxt, 
		temporalbasis* tbas ){
	rl scale=1.;
	applydeltapsi_Uinv(ncorr,nmax,deltapsi,guessxt,residxt,scale,tbas);
}


void wf::applydeltapsi_Uinv(int ncorr,int nmax, cmplx* __restrict deltapsi,
		cmplx* __restrict * __restrict guessxt,
		cmplx* __restrict * __restrict residxt,rl scale, 
		temporalbasis* tbas ){
	cmplx II=cmplx(0.,1.);
	int ot=tbas->order;
	rl dt=tbas->tmax-tbas->tmin;

	//cout << "dt in applydeltapsi_Uinv\t"<<dt<<"\n";

	//apply deltapsi vector to entry n of guessxt, adjust residuals by
	//L.Uinv.deltapsi

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

	//find H0.deltapsi, nabla.deltapsi, Olap.deltapsi
	cmplx* __restrict H0deltapsi=gbas->Hpsi(deltapsi);
	cmplx* __restrict nabladeltapsi=gbas->nablapsi(deltapsi);
	cmplx* __restrict transposenabladeltapsi=gbas->transposenablapsi(deltapsi);
	cmplx* __restrict Olapdeltapsi=gbas->Olappsi(deltapsi);
	cmplx* __restrict nablaolapdeltapsi=gbas->nablapsi(Olapdeltapsi);
	cmplx* __restrict olapnabladeltapsi=gbas->Olappsi(nabladeltapsi);
	cmplx* __restrict Olapsqdeltapsi=gbas->Olappsi(Olapdeltapsi);

	//Uinv is inverse of tridiagonal overlap matrix (less factor of dt)
	cmplx* __restrict Uinv=triD_Uinv(ot);


	//find Uinv.Qmat
	cmplx alpha=1.;
	cmplx beta=0.;
	cmplx* __restrict UinvQ=initializearray(ot*ot);
	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
			Uinv,ot,tbas->qmat,ot,&beta,UinvQ,ot);

	//find Uinv.Umat(should be diagonal w/factor dt on diagonals)
	cmplx* __restrict UinvU=initializearray(ot*ot);
	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
			Uinv,ot,tbas->umat,ot,&beta,UinvU,ot);

	//find Uinv.Amat
	cmplx* __restrict UinvA=initializearray(ot*ot);

	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
			Uinv,ot,tbas->Amat,ot,&beta,UinvA,ot);
	//find Uinv.Asqmat
	cmplx* __restrict UinvAsq=initializearray(ot*ot);

	cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ot,ot,ot,&alpha,
			Uinv,ot,tbas->Asqmat,ot,&beta,UinvAsq,ot);


	//correct residuals according to L.deltapsi
	//for tridiagonal temporal basis, umat will have bandwith of at most 2
	//q will have nonzero elements only for n<ncorr
	//Amat and Asqmat may have nonzero elements of any order
	//umat has bandwidth 2
	//for(int n=0;n<ot;n++){//)for(int n=max(0,ncorr-2);n<min(nmax,ncorr+2);n++){
	int n=ncorr;
	cmplx uval=UinvU[arrayindx(n,ncorr,ot,ot)];
	//cout << "applydeltapsi uval "<<ncorr<<" "<<n<<" "<<uval<<"\n";
	cmplx coeff=-uval;//*pow(scale,ncorr-n);
	cout << "UinvU coeff\t"<<n<<"\t"<<ncorr<<"\t"<<coeff<<"\n";
	cblas_zaxpy(nbf,&coeff,H0deltapsi,inc,residxt[n],inc);
	//}
	//qmat elts are nonzero only for n<=ncorr
	for(int n=0;n<ot;n++){//for(int n=0;n<=max(ncorr,1);n++){
		cmplx qval=UinvQ[arrayindx(n,ncorr,ot,ot)];
		//cout << "applydeltapsi qval "<<ncorr<<" "<<n<<" "<<qval<<"\n";
		cmplx coeff=II*qval;//*pow(scale,ncorr-n);
		cout << "UinvQ coeff\t"<<n<<"\t"<<ncorr<<"\t"<<coeff<<"\n";
		cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
	}
	for(int n=0;n<nmax;n++){
		cmplx Aval=UinvA[arrayindx(n,ncorr,ot,ot)];
		cmplx coeff=-II*Aval*pow(scale,ncorr-n);
		//cout << "UinvA coeff\t"<<n<<"\t"<<ncorr<<"\t"<<coeff<<"\n";
		cblas_zaxpy(nbf,&coeff,nabladeltapsi,inc,residxt[n],inc);
		cmplx Asqval=UinvAsq[arrayindx(n,ncorr,ot,ot)];
		//coeff=-AsqUval;
		coeff=Asqval*pow(scale,ncorr-n);
		//cout << "UinvAsq coeff\t"<<n<<"\t"<<ncorr<<"\t"<<coeff<<"\n";
		cblas_zaxpy(nbf,&coeff,Olapdeltapsi,inc,residxt[n],inc);
	}



	delete [] Uinv;
	delete [] UinvQ;
	delete [] UinvU;
	delete [] UinvA;
	delete [] UinvAsq;

	//delete [] AUmat;
	//delete [] AsqUmat;
	delete [] H0deltapsi;
	delete [] Olapdeltapsi;
	delete [] Olapsqdeltapsi;
	delete [] nabladeltapsi;
	delete [] transposenabladeltapsi;
	delete [] nablaolapdeltapsi;
	delete [] olapnabladeltapsi;

	}

	cmplx* __restrict * __restrict wf::
		cancelresiduals_Uinv(cmplx* __restrict * __restrict residxt,
				temporalbasis* tbas){
			//find deltaphi necessary to cancel residxt
			cout << "inside cancelresiduals_Uinv\n";
			int ot=tbas->order;

			cmplx* __restrict * __restrict deltaphi=new cmplx* __restrict[ot];
			cmplx* __restrict * __restrict tmpresid=new cmplx* __restrict[ot];
			for(int n=0;n<ot;n++){
				deltaphi[n]=initializearray(nbf);
				tmpresid[n]=arraycopy(residxt[n],nbf);
			}

			//correct residuals at levels greater than 0 (residuals at level 0 may be
			//nonzero)
			for(int n=ot-1;n>0;n--){
				//cout << "n\t"<<n<<"\n";
				cout <<"tmpresid n before = "<<n<<"\t"<<nthcomponentsignificance(tmpresid,n)<<"\n";
				cmplx* deltaphi_n=residualcorrection_Uinv(n,n,tmpresid[n],tbas);
				cout << "found deltaphi_n\t"<<sqrt(abs(gbas->norm(deltaphi_n)))<<"\n";
				applydeltapsi_Uinv(n,ot,deltaphi_n,deltaphi,tmpresid,tbas);
				cout <<"tmpresid n = "<<n<<"\t"<<nthcomponentsignificance(tmpresid,n)<<"\n";
				delete [] deltaphi_n;
			}

			//  //for every level greater than zero, correct residuals at next lowest level
			//  for(int n=1;n<ot;n++){
			//        cout <<"tmpresid n before = "<<n<<"\t"<<nthcomponentsignificance(tmpresid,n)<<"\n";
			//    cmplx* deltaphi_n=residualcorrection_Uinv(n,n-1,tmpresid[n],tbas);
			//    cout << "found deltaphi_n\n";
			//    applydeltapsi_Uinv(n,ot,deltaphi_n,deltaphi,tmpresid,tbas);
			//
			//    cout <<"tmpresid n = "<<n<<"\t"<<nthcomponentsignificance(tmpresid,n-1)<<"\n";
			//    delete [] deltaphi_n;
			//  }

			//cout << "tmpresid after correction\n";
			//for(int n=0;n<ot;n++){
			//  cout <<"norm n = "<<n<<"\t"<<nthcomponentsignificance(tmpresid,n)<<"\n";
			//}

			for(int n=0;n<ot;n++){
				delete [] tmpresid[n];
			}
			delete [] tmpresid;

			return deltaphi;
		}


	void wf::timestep_Uinv(rl tmin, rl tmax, potential* Pot,pulse* Pls,
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

		cmplx* __restrict * __restrict Uinvtargetxt=Uinvpsi(ot,targetxt);
		//scalewfxt(Uinvtargetxt,ot,dt);
		cout << "max magnitude of targetxt before scaling\n"<<maxdelta(targetxt,ot);
		scalewfxt(targetxt,ot,dt);

		//scale targetresidual by dt^-n, so that higher order terms will have
		//magnitudes comparable to lower order terms
		//scalewfxt(targetxt,ot,dt);
		cout << "max magnitude of targetxt\n"<<maxdelta(targetxt,ot);

		cmplx* __restrict * __restrict guessxt=deltaxt_Uinv(0,ot,accuracygoal,psi,
				Uinvtargetxt,tbas_triD);


		//rl enorm=errornorm(ot-1,residxt,targetxt);
		//cout <<"error norm\t"<<enorm<<"\n";
		//for(int nt=0;nt<ot;nt++){
		//  cout << nt<<"\tresidual significance\t"<<
		//    nthcomponentsignificance(residxt,targetxt,nt)<<"\n";
		//}


		for(int nt=0;nt<ot;nt++){
			cout << nt<<"\tcomponent significance\t"<<
				nthcomponentsignificance(guessxt,nt)<<"\n";
		}



		//rl truncationerror=lastcomponentsignificance(guessxt,ot);
		rl truncationerror=nthcomponentsignificance(guessxt,ot-1);
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
			//delete [] Uguessxt[n];
			delete [] guessxt[n];
			//delete [] residxt[n];
			delete [] targetxt[n];
		}
		//delete [] Uguessxt;
		delete [] guessxt;
		//delete [] residxt;
		delete [] targetxt;
	}

	cmplx* __restrict * __restrict wf::psixtdiff(
			cmplx* __restrict * __restrict psi1, 
			cmplx* __restrict * __restrict psi2, int ot){
		cmplx* __restrict * __restrict deltapsi=new cmplx*[ot];
		cmplx coeff1=1.;
		cmplx coeff2=1.;
		for(int n=0;n<ot;n++){
			deltapsi[n]=arraydiff(coeff1,psi1[n],coeff2,psi2[n],nbf);
		}
		return deltapsi;
	}

	cmplx* __restrict * __restrict 
		wf::deltaxt_Uinv(int recursionlvl, int maxrecursionlvl, rl accuracygoal,
				cmplx * __restrict psi0,
				cmplx* __restrict * __restrict targetxt,temporalbasis* tbas){
			//given initial condition psi0 and target residual tbas, find correction
			//deltaxt which cancels targetxt

			int ot=tbas->order;
			rl dt=tbas->tmax-tbas->tmin;
			cmplx* __restrict * __restrict guessxt=new cmplx* __restrict[ot];
			cmplx* __restrict * __restrict residxt=new cmplx* __restrict[ot];
			for(int n=0;n<ot;n++){
				guessxt[n]=initializearray(nbf);
				residxt[n]=initializearray(nbf);
			}

			//cout << "qmat\n";
			//printmat(tbas->qmat,ot,ot);
			//
			//cout << "umat\n";
			//printmat(tbas->umat,ot,ot);

			rl scale=scalefactor(ot,dt);

			//apply initial guess
			initialguess_Uinv(ot,psi0,guessxt,residxt,scale,tbas);
			//Upsi0(ot,psi0,guessxt,residxt,scale,tbas);


			//difference between initial residual and target must be corrected for
			//cmplx* __restrict * __restrict diffxt=psixtdiff(residxt,targetxt,ot);//psixtdiff(targetxt,residxt,ot);
			//zeroarray(diffxt[0],nbf);




			//cout << "calling Udeltaxt\n";
			cmplx* __restrict * __restrict diffxt=psixtdiff(residxt,targetxt,ot);
			//scalewfxt(Uinvdiffxt,ot,dt);
			cout << "maxdelta(diffxt)"<<maxdelta(diffxt,ot);
			//cmplx* __restrict * __restrict deltaxt=
			//  cancelresiduals_Uinv(Uinvdiffxt,tbas);//solve for U*deltapsi
			cmplx * __restrict * __restrict deltaxt=IOUinvQsolve(ot,diffxt,scale);

			cout << "maxdelta(deltaxt) after IOUinvQsolve"<<maxdelta(deltaxt,ot);

			//for(int n=0;n<ot;n++){
			//  cmplx scale=pow(dt,n);
			//  arrayscale(deltaxt[n],nbf,scale);
			//}
			//cout <<"Calling Uinvpsi\n";
			//zeroarray(deltaxt[0],nbf);//don't change initial condition
			//cout << "calling applydeltapsi\n";
			applydeltapsi_Uinv(ot,deltaxt,guessxt,residxt,scale,tbas);
			cout << "maxdelta(residxt) after applying deltapsi"<<maxdelta(residxt,ot);
			cout << "maxdelta(targetxt) after applying deltapsi"<<maxdelta(targetxt,ot);

			//target for next level of recursion is difference between 
			cmplx* __restrict * __restrict nexttarget=psixtdiff(residxt,targetxt,ot);

			cout << "recursionlevel, maxdelta(nexttarget)\t"<<recursionlvl<<"\t"<<
				maxdelta(nexttarget,ot)<<"\n";
			if((maxdelta(deltaxt,ot)>accuracygoal*1.e-6) and 
					(recursionlvl<maxrecursionlvl)){
				cmplx* __restrict nextpsi0=initializearray(nbf);
				cmplx* __restrict * __restrict nextdeltaxt=
					deltaxt_Uinv(recursionlvl+1,maxrecursionlvl,
							accuracygoal,nextpsi0,nexttarget,tbas);
				applydeltapsi(ot,nextdeltaxt,guessxt,residxt,tbas);//apply deltapsi

				for(int n=0;n<ot;n++){
					delete [] nextdeltaxt[n];
				}
				delete [] nextdeltaxt;
				delete [] nextpsi0;

			}

			//delete arrays created in this routine
			for(int n=0;n<ot;n++){
				//delete [] diffxt[n];
				delete [] residxt[n];
				delete [] deltaxt[n];
			}
			//delete [] diffxt;
			delete [] residxt;
			delete [] deltaxt;

			return guessxt;
		}

	rl wf::maxdelta(cmplx* __restrict * __restrict deltaxt,int nmax){
		rl maxdelta=0.;
		for(int nt=1;nt<nmax;nt++){
			rl delta=nthcomponentsignificance(deltaxt,nt);
			maxdelta=max(maxdelta,delta);
			//cout << nt<<"\tsignificance\t"<<delta<<"\n";
		}
		//cout << "maxdelta\t"<<maxdelta<<"\n";
		return maxdelta;
	}

	rl wf::scalefactor(int ot, rl dt){
		//in order to make backsubstitution numerically well behaved, it will be
		//useful to scale psixt[n] ->psixt[n]/pow(scalefactor,n).  This will make
		//the coefficients of the vector used in the solution comparably sized, and
		//make the backsubstitution convergent rather than divergent.  To do this,
		//we particularly want scalefactor*UinvQ[n-1,n] to be less than 0.

		//return pow(dt,1)/(2.*ot+1.);
		return 1.;//pow(dt,1.);
	}

	void wf::scalewfxt(cmplx* __restrict *__restrict wfxt,int ot, rl dt){
		for(int n=0;n<ot;n++){
			cmplx coeff=pow(scalefactor(ot,dt),-n);
			arrayscale(wfxt[n],nbf,coeff);
		}
	}

	void wf::scalewfxt_inverse(cmplx* __restrict *__restrict wfxt,int ot, rl dt){
		for(int n=0;n<ot;n++){
			cmplx coeff=pow(scalefactor(ot,dt),n);
			arrayscale(wfxt[n],nbf,coeff);
		}
	}



	cmplx* wf::Olapinvpsi(int ncorr, cmplx* __restrict resid){
		//set up and solve linear system for Olap_ij psi_j = R_i
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

		//copy olap into lhsmat
		cmplx* __restrict lhsmat=initializearray(nbf*width);
		for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
			derivedbasis* tmpbasis=gbas->elementbases[eltnum];
			int tmporder=tmpbasis->order;
			cmplx* __restrict tmpo=tmpbasis->omat;
			cmplx alpha=cmplx(0.,1.);
			cmplx beta=0.;
			//cout << "tmpo "<<eltnum<<"\n";
			//printmat(tmpo,tmporder,tmporder);
			for(int i=0;i<tmporder;i++){
				for(int j=0;j<tmporder;j++){
					int indx1=gbas->eltfirstindices[eltnum]+i;
					int indx2=gbas->eltfirstindices[eltnum]+j;
					lhsmat[bandarrayindx(indx1,indx2,ku,kl,nbf)]+=
						tmpo[arrayindx(i,j,tmporder,tmporder)];
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

	cmplx* __restrict * __restrict 
		wf::IOUinvQsolve(int ot,cmplx* __restrict * __restrict resid, rl scale){
			//solve II*Olap*Q*Uinv*dpsi=resid

			cmplx* __restrict * __restrict OinvR=new cmplx*[ot];
			for(int n=0;n<ot;n++){
				OinvR[n]=Olapinvpsi(n,resid[n]);
				//multiply by -II
				cmplx coeff=cmplx(0.,-1.);
				arrayscale(OinvR[n],nbf,coeff);
			}
			cout << "max magnitude of OinvR\n"<<maxdelta(OinvR,ot);
			cmplx* __restrict * __restrict UOinvR=Upsi(ot,OinvR,scale);
			cout << "max magnitude of UOinvR\n"<<maxdelta(UOinvR,ot);
			cmplx* __restrict * __restrict QinvUOinvR=Qinvpsi(ot,UOinvR,scale);
			cout << "max magnitude of QinvUOinvR\n"<<maxdelta(QinvUOinvR,ot);


			//delete temporary arrays
			for(int n=0;n<ot;n++){
				delete [] OinvR[n];
				delete [] UOinvR[n];
			}
			delete [] OinvR;
			delete [] UOinvR;


			return QinvUOinvR;  
		}
