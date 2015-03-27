#include "./classes.h"
#include "./wf.h"

//The routines in this file are helpers for the routines contained in wf_iterativesolution_temporaleigenvectors.cpp

cmplx* __restrict * __restrict wf::cancel0residuals_temporaleigenvector(
		int nbf,int ot,
		cmplx* __restrict resid0,
		cmplx* __restrict * __restrict VR,
		cmplx Uval,cmplx Aval, cmplx Asqval){

	int nev=ot-1;

	//find common correction which, when added to all eigenvectors, cancels
	//level 0 residual
	cmplx* __restrict deltaeigenvector=
		commoneigenvectorcorrection_level0residuals(nbf,ot,resid0,VR,Uval,Aval,
				Asqval);
	//add this correction to all eigenvectors
	cmplx* __restrict * __restrict dpsi=new cmplx*[nev];
	for(int n=0;n<ot-1;n++){
		dpsi[n]=arraycopy(deltaeigenvector,nbf);
	}
	delete [] deltaeigenvector;
	return dpsi;
}

cmplx* __restrict * __restrict wf::canceleigenvectorresiduals(
		int nbf,int ot,cmplx* __restrict evals,
		cmplx* __restrict * __restrict padrevecs,
		cmplx* __restrict * __restrict padlevecs,
		cmplx* __restrict * __restrict targetxt,temporalbasis* tbas){
	int nev=ot-1;

	cmplx* __restrict * __restrict dpsi_ev=new cmplx*[nev];
	for(int n=0;n<nev;n++){
		cmplx* __restrict targetproj_n=temporalvector_rhstarget(ot,padlevecs[n],
				targetxt);
		//cmplx scale=-1.;
		//arrayscale(targetproj_n,nbf,scale);
		cmplx Qval=evals[n];
		// Umat, Amat, Asqmat are proportional to U_R (or close, for A & Asq)
		cmplx Uval=tbas->umat[arrayindx(1,1,ot,ot)]/triD_U(1,1);
		cmplx Aval=tbas->Amat[arrayindx(1,1,ot,ot)]/triD_U(1,1);
		cmplx Asqval=tbas->Asqmat[arrayindx(1,1,ot,ot)]/triD_U(1,1);
		cmplx* __restrict dpsi_n=deltawf_eigenvector(targetproj_n,Qval,Uval,
				Aval,Asqval);
		dpsi_ev[n]=dpsi_n;
		delete [] targetproj_n;
	}
	return dpsi_ev;
}

cmplx* __restrict * __restrict 
wf::canceleigenvectorresiduals_Olap_preconditioner(
		int nbf,int ot,cmplx* __restrict evals,
		cmplx* __restrict * __restrict padrevecs,
		cmplx* __restrict * __restrict padlevecs,
		cmplx* __restrict * __restrict targetxt,temporalbasis* tbas){
	int nev=ot-1;

	cmplx* __restrict * __restrict dpsi_ev=new cmplx*[nev];
	for(int n=0;n<nev;n++){
		cmplx* __restrict targetproj_n=temporalvector_rhstarget(ot,padlevecs[n],
				targetxt);
		//cmplx scale=-1.;
		//arrayscale(targetproj_n,nbf,scale);
		cmplx Qval=evals[n];
		// Umat, Amat, Asqmat are proportional to U_R (or close, for A & Asq)
		cmplx Uval=tbas->umat[arrayindx(1,1,ot,ot)]/triD_U(1,1);
		cmplx Aval=tbas->Amat[arrayindx(1,1,ot,ot)]/triD_U(1,1);
		cmplx Asqval=tbas->Asqmat[arrayindx(1,1,ot,ot)]/triD_U(1,1);
		cmplx* __restrict dpsi_n=deltawf_eigenvector_Olap_preconditioner(
				targetproj_n,Qval,Uval,Aval,Asqval);
		dpsi_ev[n]=dpsi_n;
		delete [] targetproj_n;
	}
	return dpsi_ev;
}





cmplx* __restrict wf::
commoneigenvectorcorrection_level0residuals(int nbf,int ot,
		cmplx* __restrict resid0,
		cmplx* __restrict * __restrict VR,cmplx Uval, 
		cmplx Aval, cmplx Asqval){
	//want to find F_j such that 
	//\sum_{\alpha} O_ij F_j (Q v^{\alpha})_0- H_ij F_j (U v^{\alpha})_0=R_i0
	//ie, want to find common correction to all eigenvectors which will have the
	//effect of canceling residual at level 0

	cmplx II=cmplx(0.,1.);

	//Uval,Aval, Asqval are multiples of the tridiagonal U matrix.  Set up this
	//matrix
	cmplx* __restrict Utmp=new cmplx[ot*ot];
	cmplx* Qtmp=new cmplx[ot*ot];
	for(int n=0;n<ot;n++){
		for(int m=0;m<ot;m++){
			Utmp[arrayindx(n,m,ot,ot)]=triD_U(n,m);
			Qtmp[arrayindx(n,m,ot,ot)]=triD_Q(n,m);
		}
	}

	//for each eigenvector V, find Q.V and U.V
	cmplx* __restrict * __restrict QVR=new cmplx*[ot-1];
	cmplx* __restrict * __restrict UVR=new cmplx*[ot-1];
	for(int n=0;n<ot-1;n++){
		QVR[n]=initializearray(ot);
		UVR[n]=initializearray(ot);
		cmplx alpha=1.;
		cmplx beta=0.;
		int one=1;
		cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,&alpha,Qtmp,ot,VR[n],one,
				&beta,QVR[n],one);
		cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,&alpha,Utmp,ot,VR[n],one,
				&beta,UVR[n],one);
		//int one=1;
		//cblas_zaxpy(ot,Qtmp,VR[n],one,QVR[n],one);
		//UVR[n]=initializearray(ot);
		//cblas_zaxpy(ot,Utmp,VR[n],one,UVR[n],one);
	}

	//set up linear system obeyed by Fi
	cmplx val=-1;
	cmplx* rhsvec=arraymultiply(resid0,nbf,val);

	int kl=gbas->xorder;
	int ku=gbas->xorder;
	int width=2*kl+ku+1;

	//copy element hmat, olap into lhsmat
	cmplx* __restrict lhsmat=initializearray(nbf*width);


	for(int vecnum=0;vecnum<ot-1;vecnum++){
		cmplx vecQval=QVR[vecnum][0];
		cmplx vecUval=UVR[vecnum][0];


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

			//instead of filling in lhsmat item by item, do it with blas calls
			for(int j=0;j<tmporder;j++){
				cmplx* __restrict tmpsum=initializearray(tmporder);
				//cmplx* __restrict tmpsum2=initializearray(tmporder);

				int inc=1;
				cmplx Qcoeff=II*vecQval;
				//for arrayindx, i is the fast index, so looping through i can be
				//replaced with blas calls
				cblas_zaxpy(tmporder,&Qcoeff,tmpo+arrayindx(0,j,tmporder,tmporder),
						inc,tmpsum,inc);
				cmplx Ucoeff1=.5*Uval*vecUval;
				cblas_zaxpy(tmporder,&Ucoeff1,tmpnablasq+arrayindx(0,j,tmporder,tmporder),
						inc,tmpsum,inc);
				cmplx Ucoeff2=-Uval*vecUval;
				cblas_zaxpy(tmporder,&Ucoeff2,tmpV+arrayindx(0,j,tmporder,tmporder),
						inc,tmpsum,inc);
				cmplx Acoeff=-II*Aval*vecUval;
				cblas_zaxpy(tmporder,&Acoeff,tmpnabla+arrayindx(0,j,tmporder,tmporder),
						inc,tmpsum,inc);
				cmplx Asqcoeff=-.5*Asqval*vecUval;
				cblas_zaxpy(tmporder,&Asqcoeff,tmpo+arrayindx(0,j,tmporder,tmporder),
						inc,tmpsum,inc);


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

		}
	}


	//solve for deltapsin using zgbsv_
	int info=0;
	int neq=nbf;
	int* ipiv=new int[neq];
	int nrhs=1;
	int ldab=2*kl+ku+1;

	//cout << "calling zgbsv_\n";
	zgbsv_(&neq, &kl, &ku, &nrhs, reinterpret_cast <Lcmplx*> (lhsmat), &ldab,
			ipiv, reinterpret_cast <Lcmplx*> (rhsvec), &neq, &info);
	if(info!=0){
		cout << "problem with zgbsv_ in cancel0residuals! "<<info<<"\n";
	}

	//cout << "beginning deletions\n";
	delete [] lhsmat;
	delete [] ipiv;

	delete [] Utmp;
	delete [] Qtmp;
	for(int n=0;n<ot-1;n++){
		delete [] QVR[n];
		delete [] UVR[n];
	}
	delete [] QVR;
	delete [] UVR;

	//cout << "returning from cancel0residuals\n";
	return rhsvec;
}

cmplx* __restrict * __restrict wf::
eigenvectortowf(int nbf,int ot,
		cmplx* __restrict * __restrict dpsi,
		cmplx* __restrict * __restrict padrevecs){
	//dpsi is array of spatial vectors which are coefficients of temporal
	//eigenvectors.  This routine multiplies them out & returns them in the
	//original temporal basis.
	cmplx* __restrict * __restrict deltapsixt=new cmplx*[ot];
	for(int n=0;n<ot;n++){
		deltapsixt[n]=initializearray(nbf);
	}

	int nev=ot-1;
	for(int n=0;n<nev;n++){
		for(int m=0;m<ot;m++){
			cblas_zaxpy(nbf,&padrevecs[n][m],dpsi[n],1,deltapsixt[m],1);
		}
	}

	return deltapsixt;
}

void wf::
applydeltapsi_temporaleigenvectors(int nbf,int nmax,
		cmplx* __restrict * __restrict dpsi_ev,
		cmplx* __restrict * __restrict padrevecs,
		cmplx* __restrict * __restrict guessxt,
		cmplx* __restrict * __restrict residxt,
		temporalbasis* tbas ){
	//convert dpsi_ev from eigenvector basis to deltapsixt in normal basis, then
	//apply corrections from deltapsixt to guessxt & update residual
	cmplx * __restrict * __restrict deltapsixt=
		eigenvectortowf(nbf,nmax,dpsi_ev,padrevecs);
	applydeltapsi(nmax,deltapsixt,guessxt,residxt,tbas);

	for(int n=0;n<nmax;n++){
		delete [] deltapsixt[n];
	}
	delete [] deltapsixt;


}


void wf::
cancelresiduals_generalized_eigenvector(int nbf,int ot, cmplx* M1, cmplx* M2,
		cmplx* __restrict * __restrict guessxt,
		cmplx* __restrict * __restrict residxt,
		cmplx* __restrict * __restrict targetxt,
		temporalbasis* tbas){
	//Cancel residuals of level n>=1 using generalized eigenvector decomposition

	//find generalized eigensystem for M1[1:,1:] & M2[1:,1:]
	//test: use old QUeigensystem
	//cmplx** eigensystem=triD_QUeigensystem(ot);

	cmplx** eigensystem=generalized_eigensystem_reduced(ot,M1,M2);
	cmplx* __restrict evals=eigensystem[0];
	cmplx* __restrict levecs=eigensystem[1];
	cmplx* __restrict revecs=eigensystem[2];
	cmplx* __restrict * __restrict padrevecs=padeigenvectors(ot,revecs);
	cmplx* __restrict * __restrict padlevecs=padeigenvectors(ot,levecs);

	//cancel targetxt residuals for temporal index greater than 0
	cmplx* __restrict * __restrict dpsi_ev=
		canceleigenvectorresiduals(nbf,ot,evals,padrevecs,padlevecs,targetxt,tbas);
	//convert deltawf_eigenvectors back to tridiagonal basis & apply correction
	//to guessxt, residxt
	applydeltapsi_temporaleigenvectors(nbf,ot,dpsi_ev,padrevecs,
			guessxt,residxt,tbas);
	//delete temporary vectors
	int nev=ot-1;
	for(int n=0;n<nev;n++){
		delete [] dpsi_ev[n];
		delete [] padrevecs[n];
		delete [] padlevecs[n];
	}
	delete [] dpsi_ev;
	delete [] padrevecs;
	delete [] padlevecs;
	delete [] revecs;
	delete [] levecs;
	delete [] evals;
	delete [] eigensystem;
}

void wf::
cancelresiduals_generalized_eigenvector_Olap_preconditioner(int nbf,int ot, cmplx* M1, cmplx* M2,
		cmplx* __restrict * __restrict guessxt,
		cmplx* __restrict * __restrict residxt,
		cmplx* __restrict * __restrict targetxt,
		temporalbasis* tbas){
	//Cancel residuals of level n>=1 using generalized eigenvector decomposition

	//find generalized eigensystem for M1[1:,1:] & M2[1:,1:]
	//test: use old QUeigensystem
	//cmplx** eigensystem=triD_QUeigensystem(ot);

	//This is very similar to cancelresiduals_generalized_eigenvector, but
	//uses the overlap matrix as a right preconditioner.  Ie, instead of solving
	//L.v=R, it solves L.Oinv.O.v=R in two steps:
	//1) Solve L.w=O.R
	//2) Solve O.v=w
	//The advantage to using a spatial preconditioner is that you don't get into
	//closed cycles of convergence where one relaxation step exactly cancels the
	//other.

	int nev=ot-1;
	cmplx** eigensystem=generalized_eigensystem_reduced(ot,M1,M2);
	cmplx* __restrict evals=eigensystem[0];
	cmplx* __restrict levecs=eigensystem[1];
	cmplx* __restrict revecs=eigensystem[2];
	cmplx* __restrict * __restrict padrevecs=padeigenvectors(ot,revecs);
	cmplx* __restrict * __restrict padlevecs=padeigenvectors(ot,levecs);

	//find product of overlap matrix with target
	cmplx* __restrict * __restrict Olaptarget=new cmplx*[ot];
	for(int n=0;n<ot;n++){
		Olaptarget[n]=gbas->Olappsi(targetxt[n]);
	}



	//cancel targetxt residuals for temporal index greater than 0
	cmplx* __restrict * __restrict dpsi_ev=
		canceleigenvectorresiduals_Olap_preconditioner(nbf,ot,evals,
				padrevecs,padlevecs,
				Olaptarget,tbas);





	//convert deltawf_eigenvectors back to tridiagonal basis & apply correction
	//to guessxt, residxt
	applydeltapsi_temporaleigenvectors(nbf,ot,dpsi_ev,padrevecs,
			guessxt,residxt,tbas);
	//delete temporary vectors
	for(int n=0;n<nev;n++){
		delete [] dpsi_ev[n];
		//delete [] Olapdpsi_ev[n];
		delete [] padrevecs[n];
		delete [] padlevecs[n];
	}
	for(int n=0;n<ot;n++){
		delete [] Olaptarget[n];
	}
	delete [] dpsi_ev;
	//delete [] Olapdpsi_ev;
	delete [] Olaptarget;
	delete [] padrevecs;
	delete [] padlevecs;
	delete [] revecs;
	delete [] levecs;
	delete [] evals;
	delete [] eigensystem;
}



cmplx* __restrict * __restrict
wf::Upsi(int nmax,cmplx* __restrict * __restrict dphi){
	//find U.psi, where U does not include the factor of dt
	//for efficiency's sake, want to set this up as a series of blas calls

	//initialize and zero retmat
	cmplx* __restrict * __restrict retmat=new cmplx*[nmax];
	for(int n=0;n<nmax;n++){
		retmat[n]=initializearray(nbf);
	}


	//U is tridiagonal, except for elements U[0,2] and U[2,0];
	int one=1;
	if(nmax>=2){
		cmplx coeff=triD_U(0,2);
		cblas_zaxpy(nbf,&coeff,dphi[0],one,retmat[2],one);
		cblas_zaxpy(nbf,&coeff,dphi[2],one,retmat[0],one);
	}
	//all of the other nondiagonal elements of Umat are tridiagonal.

	//diagonal terms
	for(int n=0;n<nmax;n++){
		cmplx coeff=triD_U(n,n);
		cblas_zaxpy(nbf,&coeff,dphi[n],one,retmat[n],one);
	}

	//off diagonal terms
	for(int n=1;n<nmax;n++){
		int m=n-1;
		cmplx coeff=triD_U(n,m);
		cblas_zaxpy(nbf,&coeff,dphi[n],one,retmat[m],one);
		cblas_zaxpy(nbf,&coeff,dphi[m],one,retmat[n],one);
	}
	return retmat;
}


cmplx* __restrict * __restrict
wf::psi_swaporders(int nmax,int order1,int order2,
		cmplx* __restrict * __restrict dphi){
	//initialize and zero retmat
	cmplx* __restrict * __restrict retmat=new cmplx*[nmax];
	for(int n=0;n<nmax;n++){
		retmat[n]=arraycopy(dphi[n],nbf);
	}
	delete [] retmat[order1];
	retmat[order1]=arraycopy(dphi[order2],nbf);
	delete [] retmat[order2];
	retmat[order2]=arraycopy(dphi[order1],nbf);

	return retmat;
}


cmplx wf::residualaction(int nbf,int ot,cmplx* __restrict * __restrict guessxt,
		cmplx* __restrict * __restrict residxt){
	cmplx sum=0.;
	for(int n=0;n<ot;n++){
		int one=1;
		cmplx tmpsum=0.;
		cblas_zdotc_sub(nbf,guessxt[n],one,residxt[n],one,&tmpsum);
		sum+=tmpsum;

	}
	return sum;
}

/////////////////////////////////////////////////////////////
cmplx* __restrict wf::deltawf_eigenvector_Olap_preconditioner
(cmplx* __restrict resid, 
 cmplx Qval,cmplx Uval,
 cmplx Aval,cmplx Asqval){
	//solve Olap.(I Olap alphaQ Q - H alphaU -I p alphaA-Olap AlphaAsq ) dpsi = Olap.resid

	cmplx II=cmplx(0.,1.);
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

	//copy element hmat, olap into Lmat
	cmplx* __restrict Lmat=initializearray(nbf*width);
	for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
		derivedbasis* tmpbasis=gbas->elementbases[eltnum];
		int tmporder=tmpbasis->order;
		cmplx* __restrict tmph=tmpbasis->Hmat;
		cmplx* __restrict tmpnablasq=tmpbasis->nablasqmat;
		cmplx* __restrict tmpV=tmpbasis->Vmat;
		cmplx* __restrict tmpo=tmpbasis->omat;
		cmplx* __restrict tmpnabla=tmpbasis->nablamat;

		//instead of filling in Lmat item by item, do it with blas calls
		for(int j=0;j<tmporder;j++){
			cmplx* __restrict tmpsum=initializearray(tmporder);
			//cmplx* __restrict tmpsum2=initializearray(tmporder);

			int inc=1;
			cmplx Qcoeff=II*Qval;
			//for arrayindx, i is the fast index, so looping through i can be
			//replaced with blas calls
			cblas_zaxpy(tmporder,&Qcoeff,tmpo+arrayindx(0,j,tmporder,tmporder),
					inc,tmpsum,inc);
			cmplx Ucoeff1=.5*Uval;
			cblas_zaxpy(tmporder,&Ucoeff1,tmpnablasq+arrayindx(0,j,tmporder,tmporder),
					inc,tmpsum,inc);
			cmplx Ucoeff2=-Uval;
			cblas_zaxpy(tmporder,&Ucoeff2,tmpV+arrayindx(0,j,tmporder,tmporder),
					inc,tmpsum,inc);
			cmplx Acoeff=-II*Aval;
			cblas_zaxpy(tmporder,&Acoeff,tmpnabla+arrayindx(0,j,tmporder,tmporder),
					inc,tmpsum,inc);
			cmplx Asqcoeff=-.5*Asqval;
			cblas_zaxpy(tmporder,&Asqcoeff,tmpo+arrayindx(0,j,tmporder,tmporder),
					inc,tmpsum,inc);

			//for bandarrayindx, i is the fast index, so looping through i can be
			//replaced with blas calls
			cmplx sumcoeff=1.;
			int indx1=gbas->eltfirstindices[eltnum];
			int indx2=gbas->eltfirstindices[eltnum]+j;
			int width=2*kl+ku+1;//bandarrayindx decreases with increasing j
			cblas_zaxpy(tmporder,&sumcoeff,tmpsum,inc,
					Lmat+bandarrayindx(indx1,indx2,ku,kl,nbf),1);
			delete [] tmpsum;
		}
	}

	//set up Olap matrix in banded form
	cmplx* __restrict Omat=initializearray(nbf*width);
	for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
		derivedbasis* tmpbasis=gbas->elementbases[eltnum];
		int tmporder=tmpbasis->order;
		cmplx* __restrict tmpo=tmpbasis->omat;

		//instead of filling in Lmat item by item, do it with blas calls
		for(int j=0;j<tmporder;j++){
			cmplx coeff=1.;II*Qval;
			//for arrayindx, i is the fast index, so looping through i can be
			//replaced with blas calls
			int indx1=gbas->eltfirstindices[eltnum];
			int indx2=gbas->eltfirstindices[eltnum]+j;
			int width=2*kl+ku+1;//bandarrayindx decreases with increasing j
			cblas_zaxpy(tmporder,&coeff,tmpo+arrayindx(0,j,tmporder,tmporder),
					inc,Omat+bandarrayindx(indx1,indx2,ku,kl,nbf),inc);
		}
	}

	int kunew=ku*2;
	int klnew=kl*2;
	int widthnew=2*klnew+kunew+1;
	cmplx* OLmat=square_banded_matrix_multiply(nbf,kl,ku,Omat,kl,ku,Lmat);


	//solve for deltapsin using zgbsv_
	int info=0;
	int* ipiv=new int[neq];
	int nrhs=1;
	int ldab=2*klnew+kunew+1;

	zgbsv_(&neq, &klnew, &kunew, &nrhs, reinterpret_cast <Lcmplx*> (OLmat),
			&ldab, ipiv, reinterpret_cast <Lcmplx*> (rhsvec), &neq, &info);
	if(info!=0){
		cout << "problem with zgbsv_ in residualcorrection! "<<info<<"\n";
	}
	//delete [] AUmat;
	//delete [] AsqUmat;
	delete [] Lmat;
	delete [] Omat;
	delete [] OLmat;
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


//rl wf::lastlegendrecomponentsignificance(int ot,
//	   cmplx* __restrict * __restrict tridwf){
//  
//  cmplx coeff1=-.5;
//  cmplx coeff2=.5;
//  cmplx* __restrict lastcomponent=arraysum(coeff1,tridwf[0],
//					   coeff2,tridwf[1],nbf);
//  for(int n=2;n<ot;n++){
//    cmplx coeff=-pow(-1.,n);
//    int inc=1;
//    cblas_zaxpy(nbf,&coeff,tridwf[n],inc,lastcomponent,inc);
//  }
//  rl retval=(abs(gbas->norm(lastcomponent)))*(2./(2.*rl(ot)+1.));
//  delete [] lastcomponent;
//  return retval;
//}
rl wf::lastlegendrecomponentsignificance(int ot,rl dt,
		cmplx* __restrict * __restrict tridwf){
	cmplx* __restrict * __restrict legendrewf=triD_legendre_conversion(ot,tridwf);
	rl retval=(abs(gbas->norm(legendrewf[ot-1])))*dt*(2./(2.*(ot-1.)+1));

	for(int n=0;n<ot;n++){
		delete [] legendrewf[n];
	}
	delete legendrewf;

	return retval;
}


//cmplx* __restrict * __restrict wf
//::triD_legendre_conversion(int ot,cmplx* __restrict * __restrict tridwf){
//  cmplx* __restrict * __restrict legendrecoeffs=new cmplx*[ot];
//  //first find coeffs of P0
//  legendrecoeffs[0]=arraysum(tridwf[0],tridwf[1],nbf);
//  //find coeffs of P_{n} where cP_{n}=2*cF_{n}-cP_{n-1} (cP_{n}=coeff of P_{n})
//  for(int n=1;n<ot;n++){
//    cmplx Fcoeff=2.;
//    cmplx Pcoeff=-1.;
//    legendrecoeffs[n]=arraysum(Fcoeff,tridwf[n],Pcoeff,legendrecoeffs[n-1],nbf);
//  }
//  return legendrecoeffs;
//
//}


//cmplx* __restrict * __restrict wf
//::triD_legendre_conversion(int ot,cmplx* __restrict * __restrict tridwf){
//  //solve for legendre coefficients using tridiagonal algorithm
//
//  //copy tridwf coeffs into legendrecoeffs
//  cmplx* __restrict * __restrict darrays=new cmplx*[ot];
//  cmplx* __restrict * __restrict aarrays=new cmplx*[ot];
//  cmplx* __restrict * __restrict barrays=new cmplx*[ot];
//  cmplx* __restrict * __restrict carrays=new cmplx*[ot];
//  for(int n=0;n<ot;n++){
//    darrays[n]=arraycopy(tridwf[n],nbf);
//    aarrays[n]=initializearray(nbf);
//    barrays[n]=initializearray(nbf);
//    carrays[n]=initializearray(nbf);
//  }
//  //set up aarrays, barrays, carrays
//  
//  //aarrays are subdiagonal elements: 0.5 for n>0
//  for(int n=1;n<ot;n++){
//    for(int i=0;i<nbf;i++){
//      aarrays[n][i]=.5;
//    }
//  }
//  //barrays are diagonal elements: 0.5 for all n
//  for(int n=0;n<ot;n++){
//    for(int i=0;i<nbf;i++){
//      barrays[n][i]=.5;
//    }
//  }
//  //carrays are super diagonal elts: c_0=-.5, all others zero
//  for(int i=0;i<nbf;i++){
//      carrays[0][i]=-.5;
//  }
//
//  //tridiagonal matrix algorithm
//  //forward step: alter c_n
//  for(int i=0;i<nbf;i++){
//      carrays[0][i]/=barrays[0][i];
//  }
//  for(int n=1;n<ot;n++){
//    for(int i=0;i<nbf;i++){
//      carrays[n][i]/=(barrays[n][i]-carrays[n-1][i]*aarrays[n][i]);
//    }
//  }
//  //forward step: alter d_n
//  for(int i=0;i<nbf;i++){
//      darrays[0][i]/=barrays[0][i];
//  }
//  for(int n=1;n<ot;n++){
//    for(int i=0;i<nbf;i++){
//      darrays[n][i]=(darrays[n][i]-darrays[n-1][i]*aarrays[n][i])/
//	(barrays[n][i]-carrays[n-1][i]*aarrays[n][i]);
//    }
//  }
//
//  //backward step: backsubstitute to find coefficients
//  for(int n=ot-2;n>=0;n--){
//    for(int i=0;i<nbf;i++){
//      darrays[n][i]=darrays[n][i]-carrays[n][i]*darrays[n+1][i];
//    }
//  }
// 
//  delete [] aarrays;
//  delete [] barrays;
//  delete [] carrays;
//
//  return darrays;
//}


cmplx* __restrict * __restrict wf
::triD_legendre_conversion(int ot,cmplx* __restrict * __restrict tridwf){
	//solve for legendre coefficients using tridiagonal algorithm reverse order
	//of tridiagonal algorithm so that smaller coefficients can be solved for
	//with less error

	//copy tridwf coeffs into legendrecoeffs
	cmplx* __restrict * __restrict darrays=new cmplx*[ot];
	cmplx* __restrict * __restrict aarrays=new cmplx*[ot];
	cmplx* __restrict * __restrict barrays=new cmplx*[ot];
	cmplx* __restrict * __restrict carrays=new cmplx*[ot];
	for(int n=0;n<ot;n++){
		darrays[n]=arraycopy(tridwf[n],nbf);
		aarrays[n]=initializearray(nbf);
		barrays[n]=initializearray(nbf);
		carrays[n]=initializearray(nbf);
	}
	//set up aarrays, barrays, carrays

	//aarrays are subdiagonal elements: 0.5 for n>0
	for(int n=1;n<ot;n++){
		for(int i=0;i<nbf;i++){
			aarrays[n][i]=.5;
		}
	}
	//barrays are diagonal elements: 0.5 for all n
	for(int n=0;n<ot;n++){
		for(int i=0;i<nbf;i++){
			barrays[n][i]=.5;
		}
	}
	//carrays are super diagonal elts: c_0=-.5, all others zero
	for(int i=0;i<nbf;i++){
		carrays[0][i]=-.5;
	}

	//tridiagonal matrix algorithm
	//forward step: alter c_n
	for(int i=0;i<nbf;i++){
		aarrays[ot-1][i]/=barrays[ot-1][i];
	}
	for(int n=ot-2;n>=0;n--){
		for(int i=0;i<nbf;i++){
			aarrays[n][i]/=(barrays[n][i]-aarrays[n+1][i]*carrays[n][i]);
		}
	}
	//forward step: alter d_n
	for(int i=0;i<nbf;i++){
		darrays[ot-1][i]/=barrays[0][i];
	}
	for(int n=ot-2;n>=0;n--){
		for(int i=0;i<nbf;i++){
			darrays[n][i]=(darrays[n][i]-darrays[n+1][i]*carrays[n][i])/
				(barrays[n][i]-aarrays[n+1][i]*carrays[n][i]);
		}
	}

	//backward step: backsubstitute to find coefficients
	for(int n=ot-2;n>=0;n--){
		for(int i=0;i<nbf;i++){
			darrays[n][i]=darrays[n][i]-carrays[n][i]*darrays[n+1][i];
		}
	}

	delete [] aarrays;
	delete [] barrays;
	delete [] carrays;

	return darrays;
}





rl wf::normsum_legendre(int ot,rl dt, 
		cmplx* __restrict * __restrict legendrewf){
	rl normsum=0.;
	for(int n=0;n<ot;n++){
		rl norm=sqrt(abs(gbas->norm(legendrewf[n])));
		norm*=dt*2./(2.*n+1.);
		cout<<"norm\t"<<n<<"\t"<<norm<<"\n";
		normsum+=norm;
	}
	return normsum;
}



//cmplx* __restrict *__restrict wf
//::dipole_gauge(cmplx* __restrict *__restrict wfxt, pulse* Pls,
//		      temporalbasis* tbas){
//  //multiply wfxt by exp(II*A(t)*x) & reexpress in specetime basis
//  cmplx II=cmplx(0.,1.);
//  int ot=tbas->order;
//  int ox=gbas->xorder;
//  int ntintpts=ot;
//  int nxintpts=ox;
//
//
//  cmplx* __restrict * __restrict integralarray=new cmplx*[ot];
//  for(int n=0;n<ot;n++){
//    integralarray[n]=initializearray(nbf);
//  }
//
//  rl** tinttable=tbas->integrationtable(ntintpts);
//  rl* tintpts=tinttable[0];
//  rl* tintwts=tinttable[1];
//
//  //evaluate temporal basis functions on legendre points in t
//  cmplx* __restrict tarray=initializearray(ot*ntintpts);
//  for(int n=0;n<ntintpts;n++){
//    cmplx* tmpvec=tbas->evalbasisfuncs(tintpts[n]);
//    for(int m=0;m<ot;m++){
//      tarray[arrayindx(m,n,ot,ntintpts)]=tmpvec[m];
//    }
//    delete [] tmpvec;
//  }
//
//  cmplx* __restrict Avec=new cmplx[ntintpts];
//  cmplx* __restrict Evec=new cmplx[ntintpts];
//  for(int n=0;n<ntintpts;n++){
//    Avec[n]=Pls->Az(tintpts[n]);
//    Evec[n]=Pls->Ez(tintpts[n]);
//  }
//
//
//
//  //for every element, evaluate 
//  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
//    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
//    int eltfirstindex=gbas->eltfirstindices[eltnum];
//    //find legendre points in x
//    rl** xinttable=tmpbasis->integrationtable(nxintpts);
//    rl* xintpts=xinttable[0];
//    rl* xintwts=xinttable[1];
//    
//    //evaluate basis funcs on legendre points in x
//    cmplx* __restrict xarray=initializearray(nxintpts*ox);
//    for(int i=0;i<nxintpts;i++){
//      cmplx* tmpvec=tmpbasis->evalbasisfuncs(xintpts[i]);
//      for(int j=0;j<ox;j++){
//	xarray[arrayindx(i,j,nxintpts,ox)]=tmpvec[j];
//      }
//      delete [] tmpvec;
//    }
//
//    //read coeffs into coeffarray
//    cmplx* __restrict coeffarray=initializearray(ox*ot);
//    for(int n=0;n<ot;n++){
//      for(int i=0;i<ot;i++){
//	coeffarray[arrayindx(i,n,ox,ot)]=wfxt[n][eltfirstindex+i];
//      }
//    }
//
//
//    //find value of psi at integration points by taking the product of three
//    //matrices: psi(xi,tn)=\sum_{j,m} X_{j} (xi) . C_{jm} . T_{m}(tn)
//    //here this means xarray.coeffarray.tarray
//    //cmplx* CTarray=initializearray(ox*ntintpts);//CTarray=coeffarray.tarray
//    cmplx* XCarray=initializearray(nxintpts*ot);
//    cmplx* valarray=initializearray(nxintpts*ntintpts);
//    
//    cmplx alpha=1.;
//    cmplx beta=0.;
//
//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,nxintpts,ot,ox,
//		&alpha,xarray,nxintpts,coeffarray,ox,&beta,XCarray,nxintpts);
//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,nxintpts,ntintpts,ot,
//		&alpha,XCarray,nxintpts,tarray,ot,&beta,valarray,nxintpts);
//    
//    //now multiply each value by exp(i A(t_n) x_i)*xwt[i]*twt[n]
//    for(int i=0;i<nxintpts;i++){
//      for(int n=0;n<ntintpts;n++){
//	valarray[arrayindx(i,n,ox,ot)]*=xintpts[i]*Evec[n]*
//	  exp(-II*Avec[n]*xintpts[i])*
//	  xintwts[i]*tintwts[n];
//      }
//    }
//    delete [] XCarray;
//    
//    //now multiply (xarray_{ij})^T . valarray_{in} . (tarray_{mn})^T to get
//    //Olap.wfxt (reuse CTarray and coeffarray)
//
//    //first, find Uinv times t integrals  (multiply by t weights)
//    cmplx* __restrict tintmat=new cmplx[nxintpts*ot];
//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasTrans,nxintpts,ot,ntintpts,
//		&alpha,valarray,nxintpts,tarray,ot,&beta,tintmat,nxintpts);
////    //find Uinv.tintmat for every value of x_i
////    for(int i=0;i<nxintpts;i++){
////      cmplx* __restrict tmpvec=new cmplx[ot];
////      for(int n=0;n<ot;n++){
////	tmpvec[n]=tintmat[arrayindx(i,n,nxintpts,ot)];
////      }
////      cmplx* uinvtmpvec=tbas->oinvpsi(tmpvec);
////      for(int n=0;n<ot;n++){
////	tintmat[arrayindx(i,n,nxintpts,ot)]=uinvtmpvec[n];
////      }
////      delete [] tmpvec;
////      delete [] uinvtmpvec;
////    }
//
//    cmplx* __restrict xintmat=new cmplx[ox*ot];
//    cblas_zgemm(CblasColMajor,CblasTrans,CblasNoTrans,ox,ot,nxintpts,
//		&alpha,xarray,nxintpts,tintmat,nxintpts,&beta,xintmat,ox);
//
//    cout << "xintmat\n";
//    printmat(xintmat,ox,ot);
//
//
//    //read coeffarray into retwfxt
//    for(int i=0;i<ox;i++){
//      for(int n=0;n<ot;n++){
//	integralarray[n][eltfirstindex+i]+=xintmat[arrayindx(i,n,ox,ot)];
//      }
//    }
//
//    delete [] tintmat;
//    delete [] xintmat;
//    delete [] coeffarray;
//    delete [] valarray;
//    delete [] xintpts;
//    delete [] xintwts;
//    delete [] xinttable;
//  }
//
////  //since the previous step integrated over space and time, have to find
////  //(Oinv Uinv).integralarray before returning
////  cmplx* __restrict * __restrict retarray=new cmplx*[ot];
////
////  //find Oinv.integralarray
////  for(int n=0;n<ot;n++){
////    retarray[n]=gbas->Oinvpsi(integralarray[n]);
////  }
////
////  for(int n=0;n<ot;n++){
////    delete [] integralarray[n];
////  }
////  delete [] integralarray;
////
////  delete [] tintpts;
////  delete [] tintwts;
////  delete [] tinttable;
////
////  return retarray;
//
//  return integralarray;
//
//}

cmplx* __restrict wf::wf_gaugechange(cmplx* __restrict psi_in,
		pulse* Pls,rl t_in,rl plusminus){
	//multiply psi_in by exp(plusminus*II*A(t)*x)

	cmplx II=cmplx(0.,1.);
	int ox=gbas->xorder;
	int nxintpts=2*ox;

	rl Aval=Pls->Az(t_in);

	cmplx* integrals=initializearray(nbf);

	//for every element, evaluate 
	for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
		derivedbasis* tmpbasis=gbas->elementbases[eltnum];
		int eltfirstindex=gbas->eltfirstindices[eltnum];
		//find legendre points in x
		rl** xinttable=tmpbasis->integrationtable(nxintpts);
		rl* xintpts=xinttable[0];
		rl* xintwts=xinttable[1];

		//evaluate basis funcs on legendre points in x
		cmplx* __restrict xarray=initializearray(nxintpts*ox);
		for(int i=0;i<nxintpts;i++){
			cmplx* tmpvec=tmpbasis->evalbasisfuncs(xintpts[i]);
			for(int j=0;j<ox;j++){
				xarray[arrayindx(i,j,nxintpts,ox)]=tmpvec[j];
			}
			delete [] tmpvec;
		}

		cmplx* psivals=initializearray(nxintpts);
		cmplx alpha=1.;
		cmplx beta=0.;
		int one=1;
		cblas_zgemv(CblasColMajor,CblasNoTrans,nxintpts, ox,&alpha,xarray,nxintpts,
				psi_in+eltfirstindex,one,&beta,psivals,one);
		for(int i=0;i<nxintpts;i++){
			psivals[i]*=exp(plusminus*II*Aval*xintpts[i])*xintwts[i];
		}


		//multiply again by xarray to get integral of basis functions times
		//modified psi
		beta=1.;
		cblas_zgemv(CblasColMajor,CblasTrans,nxintpts, ox,&alpha,xarray,nxintpts,
				psivals,one,&beta,integrals+eltfirstindex,one);
		delete [] psivals;
	}

	//invert overlap matrix to get new basis function coefficients
	cmplx* retpsi=gbas->Oinvpsi(integrals);
	delete [] integrals;
	return retpsi;

}

void wf::psi_gaugechange(pulse* Pls,rl t_in,rl plusminus){
	if(psi!=0){
		cmplx* psinew=wf_gaugechange(psi,Pls,t_in,plusminus);
		delete [] psi;
		psi=psinew;
	}
}



/////////////////////////////////////////////////////////////////////////
cmplx* __restrict *__restrict wf
::dipole_gauge(cmplx* __restrict psi_in, rl En, pulse* Pls,
		temporalbasis* tbas){
	//multiply wfxt by exp(II*A(t)*x) & reexpress in specetime basis
	cmplx II=cmplx(0.,1.);
	rl lightc=137.;
	int ot=tbas->order;
	int ox=gbas->xorder;
	int ntintpts=ot*2;
	int nxintpts=ox*2;


	cmplx* __restrict * __restrict integralarray=new cmplx*[ot];
	for(int n=0;n<ot;n++){
		integralarray[n]=initializearray(nbf);
	}

	rl** tinttable=tbas->integrationtable(ntintpts);
	rl* tintpts=tinttable[0];
	rl* tintwts=tinttable[1];

	//evaluate temporal basis functions on legendre points in t
	cmplx* __restrict tarray=initializearray(ot*ntintpts);
	for(int n=0;n<ntintpts;n++){
		cmplx* tmpvec=tbas->evalbasisfuncs(tintpts[n]);
		for(int m=0;m<ot;m++){
			tarray[arrayindx(m,n,ot,ntintpts)]=tmpvec[m];
		}
		delete [] tmpvec;
	}


	rl omega=Pls->omega;
	cmplx* __restrict Avec=new cmplx[ntintpts];
	cmplx* __restrict Evec=new cmplx[ntintpts];
	for(int n=0;n<ntintpts;n++){
		Avec[n]=Pls->Az(tintpts[n]);
		Evec[n]=Pls->Ez(tintpts[n]);
	}

	//for every element, evaluate 
	for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
		//cout << "eltnum\t"<<eltnum<<"\n";
		derivedbasis* tmpbasis=gbas->elementbases[eltnum];
		int eltfirstindex=gbas->eltfirstindices[eltnum];
		//find legendre points in x
		rl** xinttable=tmpbasis->integrationtable(nxintpts);
		rl* xintpts=xinttable[0];
		rl* xintwts=xinttable[1];

		//evaluate X basis funcs on legendre points in x
		cmplx* __restrict xarray=initializearray(nxintpts*ox);
		for(int i=0;i<nxintpts;i++){
			cmplx* tmpvec=tmpbasis->evalbasisfuncs(xintpts[i]);
			for(int j=0;j<ox;j++){
				xarray[arrayindx(i,j,nxintpts,ox)]=tmpvec[j];
			}
			delete [] tmpvec;
		}
		//cout << "found xarray\n";


		cmplx* __restrict psivals=new cmplx[nxintpts];
		int inc=1;
		cmplx alpha=1.;
		cmplx beta=0.;
		cblas_zgemv(CblasColMajor,CblasNoTrans,nxintpts,ox,&alpha,xarray,nxintpts,
				psi_in+eltfirstindex,inc,&beta,psivals,inc);

		//find function values on x, t integration points & multiply by
		//gauss legendre weights
		cmplx* __restrict wfxtvals=new cmplx[nxintpts*ntintpts];
		for(int i=0;i<nxintpts;i++){
			for(int n=0;n<ntintpts;n++){
				wfxtvals[arrayindx(i,n,nxintpts,ntintpts)]=psivals[i]*
					xintpts[i]*(Evec[n])*exp(-II*Avec[n]*xintpts[i])
					*exp(-II*En*tintpts[n])
					*xintwts[i]*tintwts[n];
				//	wfxtvals[arrayindx(i,n,nxintpts,ntintpts)]+=psivals[i]*
				//	  (-pow(Avec[n],2))*exp(-II*Avec[n]*xintpts[i])
				//	  *exp(-II*En*tintpts[n])
				//	  *xintwts[i]*tintwts[n];
			}
		}
		delete [] psivals;

		//complete integration against basis functions using matrix multiplication
		//of xarray, tarray
		cmplx* __restrict Xwfxt=new cmplx[ox*ntintpts];
		cmplx* __restrict XwfxtT=new cmplx[ox*ot];
		cblas_zgemm(CblasColMajor,CblasTrans,CblasNoTrans,ox,ntintpts,nxintpts,
				&alpha,xarray,nxintpts,wfxtvals,nxintpts,&beta,Xwfxt,ox);
		cblas_zgemm(CblasColMajor,CblasNoTrans,CblasTrans,ox,ot,ntintpts,
				&alpha,Xwfxt,ox,tarray,ot,&beta,XwfxtT,ox);

		//    //test: multiply by umat
		//    cmplx* __restrict uinv=triD_Uinv(ot);
		//    cmplx* __restrict XwfxtTU=new cmplx[ox*ot];
		//    rl dt=tbas->tmax-tbas->tmin;
		//    cmplx alpha2=.25/dt;
		//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasTrans,ox,ot,ot,
		//		&alpha2,XwfxtT,ox,uinv,ot,&beta,XwfxtTU,ox);

		//read XwfxtT into integralarray
		for(int n=0;n<ot;n++){
			for(int i=0;i<ox;i++){
				integralarray[n][eltfirstindex+i]+=XwfxtT[arrayindx(i,n,ox,ot)];
			}
		}

		//    for(int n=0;n<ot;n++){
		//      cmplx* __restrict tmpOint=gbas->Oinvpsi(integralarray[n]);
		//      delete [] integralarray[n];
		//      integralarray[n]=tmpOint;
		//    }

		delete [] wfxtvals;
		delete [] xarray;
		delete [] Xwfxt;
		delete [] XwfxtT;
		//    delete [] XwfxtTU;

		delete [] xintpts;
		delete [] xintwts;
		delete [] xinttable;
	}

	delete [] Avec;
	delete [] Evec;
	delete [] tarray;
	delete [] tintpts;
	delete [] tintwts;
	delete [] tinttable;

	return integralarray;
}







///////////////////////////////////////////////////////////////////////////
//cmplx* __restrict *__restrict wf
//::timederiv_gauge(cmplx* __restrict psi_in, rl En, pulse* Pls,
//		      temporalbasis* tbas){
// //multiply wfxt by exp(II*A(t)*x) & reexpress in specetime basis
//  cmplx II=cmplx(0.,1.);
//  rl lightc=137.;
//  int ot=tbas->order;
//  int ox=gbas->xorder;
//  int ntintpts=ot*2;
//  int nxintpts=ox*2;
//
//
//  cmplx* __restrict * __restrict integralarray=new cmplx*[ot];
//  for(int n=0;n<ot;n++){
//    integralarray[n]=initializearray(nbf);
//  }
//
//  rl** tinttable=tbas->integrationtable(ntintpts);
//  rl* tintpts=tinttable[0];
//  rl* tintwts=tinttable[1];
//
// //evaluate temporal basis functions on legendre points in t
//  cmplx* __restrict tarray=initializearray(ot*ntintpts);
//  for(int n=0;n<ntintpts;n++){
//    cmplx* tmpvec=tbas->evalbasisfuncs(tintpts[n]);
//    for(int m=0;m<ot;m++){
//      tarray[arrayindx(m,n,ot,ntintpts)]=tmpvec[m];
//    }
//    delete [] tmpvec;
//  }
//
//
//  rl omega=Pls->omega;
//  cmplx* __restrict Avec=new cmplx[ntintpts];
//  cmplx* __restrict Evec=new cmplx[ntintpts];
//  for(int n=0;n<ntintpts;n++){
//    Avec[n]=Pls->Az(tintpts[n]);
//    Evec[n]=Pls->Ez(tintpts[n]);
//  }
//
//  //for every element, evaluate 
//  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
//    //cout << "eltnum\t"<<eltnum<<"\n";
//    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
//    int eltfirstindex=gbas->eltfirstindices[eltnum];
//    //find legendre points in x
//    rl** xinttable=tmpbasis->integrationtable(nxintpts);
//    rl* xintpts=xinttable[0];
//    rl* xintwts=xinttable[1];
//
//    //evaluate X basis funcs on legendre points in x
//    cmplx* __restrict xarray=initializearray(nxintpts*ox);
//    for(int i=0;i<nxintpts;i++){
//      cmplx* tmpvec=tmpbasis->evalbasisfuncs(xintpts[i]);
//      for(int j=0;j<ox;j++){
//        xarray[arrayindx(i,j,nxintpts,ox)]=tmpvec[j];
//      }
//      delete [] tmpvec;
//    }
//    //cout << "found xarray\n";
//
//
//    cmplx* __restrict psivals=new cmplx[nxintpts];
//    int inc=1;
//    cmplx alpha=1.;
//    cmplx beta=0.;
//    cblas_zgemv(CblasColMajor,CblasNoTrans,nxintpts,ox,&alpha,xarray,nxintpts,
//		psi_in+eltfirstindex,inc,&beta,psivals,inc);
//
//    //find function values on x, t integration points & multiply by
//    //gauss legendre weights
//    cmplx* __restrict wfxtvals=new cmplx[nxintpts*ntintpts];
//    for(int i=0;i<nxintpts;i++){
//      for(int n=0;n<ntintpts;n++){
//	wfxtvals[arrayindx(i,n,nxintpts,ntintpts)]=psivals[i]*
//	  ((Evec[n])*exp(-II*Avec[n]*xintpts[i]))
//	  *exp(-II*En*tintpts[n])
//	  *xintwts[i]*tintwts[n];
//      }
//    }
//    delete [] psivals;
//    
//    //complete integration against basis functions using matrix multiplication
//    //of xarray, tarray
//    cmplx* __restrict Xwfxt=new cmplx[ox*ntintpts];
//    cmplx* __restrict XwfxtT=new cmplx[ox*ot];
//    cblas_zgemm(CblasColMajor,CblasTrans,CblasNoTrans,ox,ntintpts,nxintpts,
//		&alpha,xarray,nxintpts,wfxtvals,nxintpts,&beta,Xwfxt,ox);
//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasTrans,ox,ot,ntintpts,
//		&alpha,Xwfxt,ox,tarray,ot,&beta,XwfxtT,ox);
//
//    //test: multiply by UinvQ
//    cmplx* __restrict uinv=triD_Uinv(ot);
//    cmplx* __restrict XwfxtTU=new cmplx[ox*ot];
//    cmplx* __restrict XwfxtTUQ=new cmplx[ox*ot];
//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,ox,ot,ot,
//		&alpha,XwfxtT,ox,uinv,ot,&beta,XwfxtTU,ox);
//    rl dt=tbas->tmax-tbas->tmin;
//    cmplx alpha2=alpha/dt;
//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasTrans,ox,ot,ot,
//		&alpha,XwfxtT,ox,tbas->qmat,ot,&beta,XwfxtTUQ,ox);
//    
//    //read XwfxtT into integralarray
//    for(int n=0;n<ot;n++){
//      for(int i=0;i<ox;i++){
//	integralarray[n][eltfirstindex+i]+=XwfxtTUQ[arrayindx(i,n,ox,ot)];
//      }
//    }
//
////    for(int n=0;n<ot;n++){
////      cmplx* __restrict tmpOint=gbas->Oinvpsi(integralarray[n]);
////      delete [] integralarray[n];
////      integralarray[n]=tmpOint;
////    }
//
//    delete [] wfxtvals;
//    delete [] xarray;
//    delete [] Xwfxt;
//    delete [] XwfxtT;
//    delete [] XwfxtTU;
//
//    delete [] xintpts;
//    delete [] xintwts;
//    delete [] xinttable;
//  }
//
//  delete [] Avec;
//  delete [] Evec;
//  delete [] tarray;
//  delete [] tintpts;
//  delete [] tintwts;
//  delete [] tinttable;
//
//  return integralarray;
//}

