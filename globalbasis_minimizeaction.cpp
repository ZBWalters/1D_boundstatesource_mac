#include "./classes.h"
//#include "./globalbasis.h"
#include "./linearsystemsetup.h"
//#include "./temporalbasis.h"


cmplx* globalbasis::minimizeaction(cmplx* psi_global, rl t1, rl t2){
	//set up and solve linear system corresponding to time step from time t1 to t2

	//  //first, need up update temporal matrices to correspond to time step
	//  cmplx** gls_band=globallinearsystem_banded(psi_global);
	//  cmplx* gls_rhs=gls_band[0];
	//  cmplx* gls_lhs=gls_band[1];
	//
	//  int nxbases=nelts*xorder-(nelts-1);//number of x basis functions, less
	//                                     //double counted boundary functions
	//  int neq=nxbases*(torder-1);
	//  int nrhs=1;
	//  int neqelt=xorder*(torder-1);//number of equations per element
	//  int ku=neqelt;
	//  int kl=neqelt;
	//  int ldab=2*kl+ku+1;
	//  int* ipiv=new int[neq];
	//  int info=0;
	//  zgbsv_(&neq, &kl, &ku, &nrhs,gls_lhs,&ldab,ipiv,gls_rhs,&neq,&info);
	//  if(info!=0){
	//    cout << "problem with zgbsv_ in timestep! "<<info<<"\n";
	//  }


	//set up and solve linear system corresponding to a single timestep using
	//lagrange multipliers to specify initial conditions
	cmplx** gls_band=globallinearsystem_lagrange_banded(psi_global);
	cmplx* gls_rhs=gls_band[0];
	cmplx* gls_lhs=gls_band[1];


	int nxbases=nelts*xorder-(nelts-1);//number of x basis functions, less
	//double counted boundary functions
	int neq=nxbases*(torder+1);
	int nrhs=1;
	int neqelt=xorder*(torder+1);//number of equations per element
	int ku=neqelt;
	int kl=neqelt;
	int ldab=2*kl+ku+1;
	int* ipiv=new int[neq];
	int info=0;



	zgbsv_(&neq,  &kl,  &ku,  &nrhs, reinterpret_cast <Lcmplx*> (gls_lhs), &ldab,
			ipiv, reinterpret_cast <Lcmplx*> (gls_rhs), &neq, &info);
	if(info!=0){
		cout << "problem with zgbsv_ in timestep! "<<info<<"\n";
	}

	//copy solution into wfxt
	cmplx* __restrict psixt=initializearray(nxbases*torder);
	for(int i=0;i<nxbases;i++){
		for(int n=0;n<torder;n++){
			psixt[arrayindx(i,n,nxbases,torder)]=gls_rhs[lsindx(i,n,nxbases,torder+1)];
		}
	}


	delete [] ipiv;
	delete [] gls_rhs;
	delete [] gls_lhs;
	delete [] gls_band;
	return psixt;
}


cmplx* globalbasis::minimizeactioncorrection(cmplx* psixt_global, 
		cmplx* dpsi0){
	//set up and solve linear system corresponding to time step from time t1 to
	//t2 initial guess for minimum action solution is psixt_global initial
	//condition for correction is dpsi0 returns correction to psixt_global which
	//minimizes action, subject to constraint that correction must have initial
	//condition dpsi0

	//set up and solve linear system corresponding to a single timestep using
	//lagrange multipliers to specify initial conditions
	cmplx** gls_band=globallinearsystem_lagrange_banded(psixt_global,dpsi0);
	cmplx* gls_rhs=gls_band[0];
	cmplx* gls_lhs=gls_band[1];

	int nxbases=nelts*xorder-(nelts-1);//number of x basis functions, less
	//double counted boundary functions
	int neq=nxbases*(torder+1);
	int nrhs=1;
	int neqelt=xorder*(torder+1);//number of equations per element
	int ku=neqelt;
	int kl=neqelt;
	int ldab=2*kl+ku+1;
	int* ipiv=new int[neq];
	int info=0;

	cout << "gls_rhs\n";
	printmat(gls_rhs,neq,1);

	cout << "test whether gls_rhs should be multiplied by -1\n";
	for(int i=0;i<neq;i++){
		gls_rhs[i]*=-1.;
	}

	//cout << "gls_lhs\n";
	//printmat(gls_lhs,ldab,neq);

	zgbsv_(&neq,  &kl,  &ku,  &nrhs, reinterpret_cast <Lcmplx*> (gls_lhs), &ldab,
			ipiv, reinterpret_cast <Lcmplx*> (gls_rhs), &neq, &info);
	if(info!=0){
		cout << "problem with zgbsv_ in timestep! "<<info<<"\n";
	}

	cout << "gls solution\n";
	printmat(gls_rhs,neq,1);

	//copy solution into wfxt
	cmplx* __restrict dpsixt=initializearray(nxbases*torder);
	for(int i=0;i<nxbases;i++){
		for(int n=0;n<torder;n++){
			dpsixt[arrayindx(i,n,nxbases,torder)]=gls_rhs[lsindx(i,n,nxbases,torder+1)];
		}
	}
	cout <<"dpsixt\n";
	printmat(dpsixt,nxbases,torder);


	delete [] ipiv;
	delete [] gls_rhs;
	delete [] gls_lhs;
	delete [] gls_band;
	return dpsixt;
}

cmplx* globalbasis::eltLpsi(int eltnum){
	return eltLpsi(eltnum,psi_xt);
}

cmplx* globalbasis::eltLpsi(int eltnum,cmplx* psi_xt_in){
	derivedbasis* tmpbasis=elementbases[eltnum];

	int tmpxorder=tmpbasis->order;
	int firstindx=arrayindx(eltfirstindices[eltnum],0,nbf,torder);

	cmplx** xmats=new cmplx*[3];
	xmats[0]=tmpbasis->omat;
	xmats[1]=tmpbasis->nablamat;
	xmats[2]=tmpbasis->nablasqmat;

	cmplx** tmats=new cmplx*[3];
	tmats[0]=tbasis->umat;
	tmats[1]=tbasis->qmat;
	tmats[2]=tbasis->Amat;

	cmplx* __restrict tmppsi=new cmplx[tmpxorder*torder];
	for(int i=0;i<tmpxorder;i++){
		for(int n=0;n<torder;n++){
			tmppsi[arrayindx(i,n,tmpxorder,torder)]=
				psi_xt_in[arrayindx(firstindx+i,n,nbf,torder)];
		}
	}

	cmplx* retpsi=Lpsi(tmpbasis->order,torder,xmats,tmats,
			tmppsi);

	cout << "elt tmppsi\n";
	printmat(tmppsi,tmpxorder,torder);

	cout << "elt Lpsi "<<eltnum<<"\n";
	printmat(retpsi,tmpxorder,torder);
	delete [] tmppsi;
	delete [] xmats;
	delete [] tmats;
	return retpsi;

}

cmplx* globalbasis::globalLpsi(){
	return globalLpsi(psi_xt);
}

cmplx* globalbasis::globalLpsi(cmplx* psi_xt_in){
	cmplx* psiret=initializearray(nbf*torder);
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		cmplx* tmpLpsi=eltLpsi(eltnum,psi_xt_in);
		for(int i=0;i<tmpbasis->order;i++){
			for(int n=0;n<torder;n++){
				psiret[arrayindx(eltfirstindices[eltnum]+i,n,nbf,torder)]+=
					tmpLpsi[arrayindx(i,n,tmpbasis->order,torder)];
			}
		}
		delete [] tmpLpsi;
	}
	cout << "global Lpsi\n";
	printmat(psiret,nbf,torder);
	return psiret;
}

cmplx globalbasis::psiLpsi(cmplx* psi_xt_1,cmplx* psi_xt_2){
	cmplx* Lpsi2=globalLpsi(psi_xt_2);
	cmplx* dot=new cmplx[1];
	int ncoeff=nbf*torder;
	cmplx alpha=cmplx(1.,0.);
	cmplx beta=0.;
	cblas_zgemv(CblasColMajor,CblasConjTrans,ncoeff,1,&alpha,psi_xt_1,ncoeff,Lpsi2,1,
			&beta,dot,1);
	delete [] Lpsi2;
	cmplx retval=dot[0];
	delete [] dot;
	//cout << "psiLpsi "<<retval<<"\n";
	return retval;
}

cmplx globalbasis::action(){
	return action(psi_xt);
}

cmplx globalbasis::action(cmplx* psi_xt_in){
	return psiLpsi(psi_xt_in,psi_xt_in);
}
