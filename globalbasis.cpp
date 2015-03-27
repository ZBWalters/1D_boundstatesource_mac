#include "./classes.h"
#include "./globalbasis.h"//includes derivedbasis.h

globalbasis::globalbasis(int nelts_in, int nbdyelts_in,
		rl dx,rl dt, int xorder_in,int torder_in, 
		rl laguerrekappa_in, rl ecstheta_in){
	//cout << "entering globalbasis\n";
	nelts=nelts_in;
	nbdyelts=nbdyelts_in;
	xorder=xorder_in;
	torder=torder_in;
	laguerrekappa=laguerrekappa_in;
	ecstheta=ecstheta_in;

	//cout << "setting up legendrebases\n";
	primitivetbasis=new legendrebasis(torder);
	primitivexbasis=new legendrebasis(xorder);
	//cout << "setting up laguerre bases\n";
	//cout << "laguerre kappa "<<laguerrekappa<<"\n";
	primitivebdybasis_l=new laguerrebasis(xorder,abs(laguerrekappa));
	primitivebdybasis_r=new laguerrebasis(xorder,-abs(laguerrekappa));
	//cout << "primitive boundary kappa l "<<primitivebdybasis_l->kappa<<"\n";
	//cout << "primitive boundary kappa r "<<primitivebdybasis_r->kappa<<"\n";


	//cout << "set up primitivexbasis --omat\n";
	//primitivexbasis->printomat();
	//cout << "nablamat\n";
	//primitivexbasis->printnablamat();
	//cout << "nablasqmat\n";
	//primitivexbasis->printnablasqmat();


	//For both the temporal and spatial basis, it is useful to use linear
	//combinations of the primitive bases rather than the primitive bases
	//directly.  Thus, we must set up the basis change matrices.

	//cout << "setting up boundary & initial condition bases\n";
	//cmplx* icbasis_legendre=leftborderfunctionbasis(torder);
	//cout << "icbasis\n";
	cmplx* identitybasis_legendre=identitybasis(torder);
	cmplx* bcbasis_legendre=borderfunctionbasis(xorder);
	//cout << "bcbasis\n";
	//printmat(bcbasis_legendre,xorder,xorder);
	cmplx* lbasis_laguerre=rightborderfunctionbasis_laguerre(xorder);
	//cout << "lbasis_laguerre\n";
	cmplx* rbasis_laguerre=leftborderfunctionbasis_laguerre(xorder);
	//cout << "rbasis_laguerre\n";

	//set up tbasis
	tbasis=new temporalbasis(primitivetbasis,identitybasis_legendre,0.,dt);

	//cout << "global basis constructor tmats\n";
	//cout << "umat\n";
	//printmat(tbasis->umat,torder,torder);
	//cout << "qmat\n";
	//printmat(tbasis->qmat,torder,torder);
	//cout << "Amat\n";
	//printmat(tbasis->Amat,torder,torder);



	//set up array of element boundaries
	elementbdys=new rl[nelts-1];
	for(int i=0;i<nelts-1;i++){
		elementbdys[i]=-dx*(rl(nelts)-2.)/2.+rl(i)*dx;
		//cout <<"elementbdy "<<i<<" "<<elementbdys[i]<<"\n";
	}
	//cout<<"elementbdys\n";
	//printmat(elementbdys,1,nelts-1);

	//set up array of local bases
	elementbases=new derivedbasis*[nelts];
	for(int i=0;i<nelts;i++){
		elementbases[i]=0;
	}
	//first and last element bases are elements derived from laguerre basis
	//intermediate element bases are elements derived from legendre basis
	elementbases[0]=new derivedbasis(primitivebdybasis_l,lbasis_laguerre,xorder,
			elementbdys[0]);
	elementbases[nelts-1]=new derivedbasis(primitivebdybasis_r,rbasis_laguerre,
			xorder,elementbdys[nelts-2]);


	for(int i=1;i<nelts-1;i++){
		elementbases[i]=new derivedbasis(primitivexbasis,bcbasis_legendre,xorder,
				elementbdys[i-1],elementbdys[i]);
		//cout << "set up elementbases "<<i<<"\n";
		//cout << "olap "<<i<<"\n";
		//elementbases[i]->printomat();
		//cout << "nabla "<<i<<"\n";
		//elementbases[i]->printnablamat();
		//cout << "nablasq "<<i<<"\n";
		//elementbases[i]->printnablasqmat();

	}
	//cout << "back from derived basis setup";
	//for(int i=1;i<nelts-1;i++){
	//  cout << i<< " xmin "<<elementbases[i]->xmin << "\n";
	//  cout << "olap "<<i<<"\n";
	//  elementbases[i]->printomat();
	//  cout << "\n";
	//}

	//eltfirstindices corresponds to the index of the global function basis
	//which corresponds to a particular element's zeroth basis function
	eltfirstindices=new int[nelts];
	eltfirstindices[0]=0;
	for(int i=1;i<nelts;i++){
		eltfirstindices[i]=eltfirstindices[i-1]+elementbases[i]->order-1;
	}

	//allocate psi, psi_xt
	//total number of basis functions = sum of orders in each element, minus
	//(nelts-1) to account for double counting.
	nbf=0;
	for(int i=0;i<nelts;i++){
		nbf+=elementbases[i]->order;
	}
	nbf-=(nelts-1);

	psi=new cmplx[nbf];
	psi_xt=new cmplx[nbf*torder];



	//delete [] icbasis_legendre;
	delete [] identitybasis_legendre;
	delete [] bcbasis_legendre;
	delete [] lbasis_laguerre;
	delete [] rbasis_laguerre;

}


globalbasis::globalbasis(potential* pot,pulse* Pls, rl En_kappa, rl deltaphi, 
		rl ecsbdy_in, rl ecsbuffer,rl dt, int xorder_in,
		int torder_in, 
		rl laguerrekappa_in, rl ecstheta_in){
	//cout << "entering globalbasis\n";
	//nelts=nelts_in;
	//nbdyelts=nbdyelts_in;
	xorder=xorder_in;
	torder=torder_in;
	laguerrekappa=laguerrekappa_in;
	ecstheta=ecstheta_in;
	ecsbdy=ecsbdy_in;

	int nelts_zerobased=0;
	int necselts_zerobased=0;
	//find element boundaries for semi-infinite range [0,infinity)
	rl* eltbdys_zerobased=
		elementboundaries_zerobased(pot,Pls,En_kappa,deltaphi,ecsbdy,ecsbuffer,ecstheta,
				nelts_zerobased,necselts_zerobased);

	//cout << "eltbdys_zerobased\n";
	//printmat(eltbdys_zerobased,nelts_zerobased,1);

	//copy element boundaries from [0,infinity] to fully infinite range
	//[-infinity, infinity]
	nelts=2*(nelts_zerobased-1);
	nbdyelts=necselts_zerobased+1;//+1 accounts for end cap elements
	elementbdys=new rl[nelts+1];
	int zeroindx=nelts_zerobased-1;
	for(int i=0;i<nelts_zerobased;i++){
		int indx1=zeroindx+i;
		elementbdys[indx1]=eltbdys_zerobased[i];
		int indx2=zeroindx-i;
		elementbdys[indx2]=-eltbdys_zerobased[i];
	}
	//cout << "elementbdys\n";
	//printmat(elementbdys,nelts+1,1);
	//cout << "after elementbdys\n";
	delete [] eltbdys_zerobased;

	//cout << "setting up legendrebases\n";
	primitivetbasis=new legendrebasis(torder);
	primitivexbasis=new legendrebasis(xorder);
	//cout << "setting up laguerre bases\n";
	//cout << "laguerre kappa "<<laguerrekappa<<"\n";
	primitivebdybasis_l=new laguerrebasis(xorder,abs(laguerrekappa));
	primitivebdybasis_r=new laguerrebasis(xorder,-abs(laguerrekappa));
	//cout << "primitive boundary kappa l "<<primitivebdybasis_l->kappa<<"\n";
	//cout << "primitive boundary kappa r "<<primitivebdybasis_r->kappa<<"\n";


	//cout << "set up primitivexbasis --omat\n";
	//primitivexbasis->printomat();
	//cout << "nablamat\n";
	//primitivexbasis->printnablamat();
	//cout << "nablasqmat\n";
	//primitivexbasis->printnablasqmat();


	//For both the temporal and spatial basis, it is useful to use linear
	//combinations of the primitive bases rather than the primitive bases
	//directly.  Thus, we must set up the basis change matrices.

	//cout << "setting up boundary & initial condition bases\n";
	//cmplx* icbasis_legendre=leftborderfunctionbasis(torder);
	//cout << "icbasis\n";
	cmplx* identitybasis_legendre=identitybasis(torder);
	cmplx* bcbasis_legendre=borderfunctionbasis(xorder);
	//cout << "bcbasis\n";
	//printmat(bcbasis_legendre,xorder,xorder);
	cmplx* lbasis_laguerre=rightborderfunctionbasis_laguerre(xorder);
	//cout << "lbasis_laguerre\n";
	cmplx* rbasis_laguerre=leftborderfunctionbasis_laguerre(xorder);
	//cout << "rbasis_laguerre\n";

	//set up tbasis
	tbasis=new temporalbasis(primitivetbasis,identitybasis_legendre,0.,dt);

	//cout << "global basis constructor tmats\n";
	//cout << "umat\n";
	//printmat(tbasis->umat,torder,torder);
	//cout << "qmat\n";
	//printmat(tbasis->qmat,torder,torder);
	//cout << "Amat\n";
	//printmat(tbasis->Amat,torder,torder);



	//  //set up array of element boundaries
	//  elementbdys=new rl[nelts-1];
	//  for(int i=0;i<nelts-1;i++){
	//    elementbdys[i]=-dx*(rl(nelts)-2.)/2.+rl(i)*dx;
	//    //cout <<"elementbdy "<<i<<" "<<elementbdys[i]<<"\n";
	//  }
	//  //cout<<"elementbdys\n";
	//  //printmat(elementbdys,1,nelts-1);

	//set up array of local bases
	elementbases=new derivedbasis*[nelts];
	for(int i=0;i<nelts;i++){
		elementbases[i]=0;
	}
	//first and last element bases are elements derived from laguerre basis
	//intermediate element bases are elements derived from legendre basis
	elementbases[0]=new derivedbasis(primitivebdybasis_l,lbasis_laguerre,xorder,
			elementbdys[0]);
	elementbases[nelts-1]=new derivedbasis(primitivebdybasis_r,rbasis_laguerre,
			xorder,elementbdys[nelts]);//elementbdys[nelts-2]);


	for(int i=1;i<nelts-1;i++){
		elementbases[i]=new derivedbasis(primitivexbasis,bcbasis_legendre,xorder,
				elementbdys[i-1],elementbdys[i]);
		//cout << "set up elementbases "<<i<<"\n";
		//cout << "olap "<<i<<"\n";
		//elementbases[i]->printomat();
		//cout << "nabla "<<i<<"\n";
		//elementbases[i]->printnablamat();
		//cout << "nablasq "<<i<<"\n";
		//elementbases[i]->printnablasqmat();

	}

	cout << "ecs nbdyelts\t"<<nbdyelts<<"\n";
	//ecssetup(nbdyelts,ecstheta);
	ecssetup(nbdyelts,ecstheta,pot);


	//cout << "back from derived basis setup";
	//for(int i=1;i<nelts-1;i++){
	//  cout << i<< " xmin "<<elementbases[i]->xmin << "\n";
	//  cout << "olap "<<i<<"\n";
	//  elementbases[i]->printomat();
	//  cout << "\n";
	//}

	//eltfirstindices corresponds to the index of the global function basis
	//which corresponds to a particular element's zeroth basis function
	eltfirstindices=new int[nelts];
	eltfirstindices[0]=0;
	for(int i=1;i<nelts;i++){
		eltfirstindices[i]=eltfirstindices[i-1]+elementbases[i]->order-1;
	}

	//allocate psi, psi_xt
	//total number of basis functions = sum of orders in each element, minus
	//(nelts-1) to account for double counting.
	nbf=0;
	for(int i=0;i<nelts;i++){
		nbf+=elementbases[i]->order;
	}
	nbf-=(nelts-1);

	psi=new cmplx[nbf];
	psi_xt=new cmplx[nbf*torder];



	//delete [] icbasis_legendre;
	delete [] identitybasis_legendre;
	delete [] bcbasis_legendre;
	delete [] lbasis_laguerre;
	delete [] rbasis_laguerre;

}



globalbasis::~globalbasis(){
	for(int i=0;i<nelts;i++){
		delete elementbases[i];
	}
	//delete elementbases[0];
	//delete elementbases[1];
	//delete elementbases[2];
	delete [] elementbases;
	delete [] elementbdys;
	delete primitivexbasis;
	delete primitivetbasis;
	delete primitivebdybasis_l;
	delete primitivebdybasis_r;
	delete tbasis;
	delete []  eltfirstindices;
	delete [] psi;
	delete [] psi_xt;


}

void globalbasis::setupVmats(potential* pot){
	//cout << "setupVmats nelts "<<nelts<<"\n";
	for(int i=0;i<nelts;i++){
		//cout << "element bases "<<i<<"\n";
		elementbases[i]->Vmatsetup(pot);
		//cout << "Vmat\n";
		//printmat(elementbases[i]->Vmat,elementbases[i]->order,elementbases[i]->order);
	}
}

rl globalbasis::norm(cmplx* psi_in){
	return innerproduct(psi_in,psi_in);
}

rl globalbasis::innerproduct(cmplx* psi1, cmplx* psi2){
	rl sum=0.;
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		int tmpindx=eltfirstindices[eltnum];
		rl tmpsum=tmpbasis->innerproduct(&psi1[tmpindx],&psi2[tmpindx]);
		//cout << "tmpsum\t"<<tmpsum<<"\n";
		sum+=tmpsum;
	}
	//cout << "sum norm\t"<<sum<<"\n";
	return sum;

}

//rl globalbasis::innerproduct(cmplx* psi1, cmplx* psi2){
//  cmplx* psitmp=Olappsi(psi2);
//  int one=1;
//  cmplx retval=0.;
//  cblas_zdotc_sub(nbf,psi1,one,psitmp,one,&retval);
//  delete [] psitmp;
//  return abs(retval);
//}

cmplx* globalbasis::Olappsi(cmplx* psi_in){
	//returns Omat.psi for the given wf
	cmplx* __restrict psiret=initializearray(nbf);//new cmplx[nbf];
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		int indx1=eltfirstindices[eltnum];
		cmplx* __restrict psitmp=tmpbasis->Olappsi(&psi_in[indx1]);
		for(int i=0;i<tmpbasis->order;i++){
			psiret[indx1+i]+=psitmp[i];
		}
		delete [] psitmp;
	}
	return psiret;
}

cmplx* globalbasis::Oinvpsi(cmplx* psi_in){
	cmplx* rhsvec=arraycopy(psi_in,nbf);

	int kl=xorder;
	int ku=xorder;
	int width=2*kl+ku+1;

	cmplx* __restrict lhsmat=initializearray(nbf*width);
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		int tmporder=tmpbasis->order;
		cmplx* __restrict tmpo=tmpbasis->omat;
		cmplx coeff=1.;
		for(int j=0;j<tmporder;j++){
			int indx1=eltfirstindices[eltnum];
			int indx2=eltfirstindices[eltnum]+j;
			int inc=1;
			cblas_zaxpy(tmporder,&coeff,tmpo+arrayindx(0,j,tmporder,tmporder),
					inc,lhsmat+bandarrayindx(indx1,indx2,ku,kl,nbf),inc);
		}
	}
	//solve for Oinvpsi using zgbsv_
	int info=0;
	int* ipiv=new int[nbf];
	int nrhs=1;
	int ldab=2*kl+ku+1;
	zgbsv_(&nbf, &kl, &ku, &nrhs, reinterpret_cast <Lcmplx*> (lhsmat), &ldab,
			ipiv, reinterpret_cast <Lcmplx*> (rhsvec), &nbf, &info);
	if(info!=0){
		cout << "problem with zgbsv_ in residualcorrection! "<<info<<"\n";
	}
	delete [] lhsmat;
	delete [] ipiv;

	return rhsvec;
}

cmplx* globalbasis::Hpsi(cmplx* psi_in){
	//returns Hmat.psi for the given wf by calculating H.psi in each element &
	//adding them together
	cmplx* __restrict psiret=initializearray(nbf);//new cmplx[nbf];
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		int indx1=eltfirstindices[eltnum];
		cmplx* __restrict psitmp=tmpbasis->Hpsi(&psi_in[indx1]);
		for(int i=0;i<tmpbasis->order;i++){
			psiret[indx1+i]+=psitmp[i];
		}
		delete [] psitmp;
	}
	return psiret;
}

cmplx* globalbasis::Vpsi(cmplx* psi_in){
	//returns Hmat.psi for the given wf by calculating H.psi in each element &
	//adding them together
	cmplx* __restrict psiret=initializearray(nbf);//new cmplx[nbf];
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		int indx1=eltfirstindices[eltnum];
		cmplx* __restrict psitmp=tmpbasis->Vpsi(&psi_in[indx1]);
		for(int i=0;i<tmpbasis->order;i++){
			psiret[indx1+i]+=psitmp[i];
		}
		delete [] psitmp;
	}
	return psiret;
}

cmplx* globalbasis::ppsi(cmplx* psi_in){
	cmplx II=cmplx(0.,1.);
	cmplx* psiret=nablapsi(psi_in);
	arrayscale(psiret,nbf,II);
	return psiret;
}

cmplx* globalbasis::nablapsi(cmplx* psi_in){
	//returns nablamat.psi for the given wf by calculating H.psi in each element
	//& adding them together
	cmplx* __restrict psiret=initializearray(nbf);//new cmplx[nbf];
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		int indx1=eltfirstindices[eltnum];
		cmplx* __restrict psitmp=tmpbasis->Xmatpsi(tmpbasis->nablamat,&psi_in[indx1]);
		for(int i=0;i<tmpbasis->order;i++){
			psiret[indx1+i]+=psitmp[i];
		}
		delete [] psitmp;
	}
	return psiret;
}

cmplx* globalbasis::nablasqpsi(cmplx* psi_in){
	//returns nablamat.psi for the given wf by calculating H.psi in each element
	//& adding them together
	cmplx* __restrict psiret=initializearray(nbf);//new cmplx[nbf];
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		int indx1=eltfirstindices[eltnum];
		cmplx* __restrict psitmp=tmpbasis->Xmatpsi(tmpbasis->nablasqmat,&psi_in[indx1]);
		for(int i=0;i<tmpbasis->order;i++){
			psiret[indx1+i]+=psitmp[i];
		}
		delete [] psitmp;
	}
	return psiret;
}

cmplx* globalbasis::transposenablapsi(cmplx* psi_in){
	//returns nablamat.psi for the given wf by calculating H.psi in each element
	//& adding them together
	cmplx* __restrict psiret=initializearray(nbf);//new cmplx[nbf];
	for(int eltnum=0;eltnum<nelts;eltnum++){
		derivedbasis* tmpbasis=elementbases[eltnum];
		int indx1=eltfirstindices[eltnum];
		int tmporder=tmpbasis->order;
		cmplx* __restrict nablatrans=transpose(tmpbasis->nablamat,tmporder,tmporder);
		cmplx* __restrict psitmp=tmpbasis->Xmatpsi(nablatrans,&psi_in[indx1]);
		delete [] nablatrans;
		for(int i=0;i<tmpbasis->order;i++){
			psiret[indx1+i]+=psitmp[i];
		}
		delete [] psitmp;
	}
	return psiret;
}


void globalbasis::updateHmats(potential* pot){
	for(int eltnum=0;eltnum<nelts;eltnum++){
		elementbases[eltnum]->updateHmat(pot);
	}
}



void globalbasis::ecssetup(int nbdyelts_in, rl ecstheta_in){
	nbdyelts=nbdyelts_in;
	ecstheta=ecstheta_in;

	for(int i=0;i<nbdyelts;i++){
		//cout << "rotating ecs elt\t"<<i<<"\n";
		elementbases[i]->ecsrotate(ecstheta);
		elementbases[nelts-1-i]->ecsrotate(ecstheta);
	}
}

void globalbasis::ecssetup(int nbdyelts_in, rl ecstheta_in,potential* pot){
	nbdyelts=nbdyelts_in;
	ecstheta=ecstheta_in;

	for(int i=0;i<nbdyelts;i++){
		//cout << "rotating ecs elt\t"<<i<<"\n";
		elementbases[i]->ecsrotate(ecstheta,pot,ecsbdy,ecstheta);
		elementbases[nelts-1-i]->ecsrotate(ecstheta,pot,ecsbdy,ecstheta);
	}
}

int globalbasis::eltnum(rl x){
	//return the index of the element corresponding to a particular x value
	//use a simple binary search to find the correct index
	if(elementbdys==0){
		cout << "null elementbdys!\n";
	}
	//cout << "inside eltnum\t"<<x<<"\n";
	//cout <<"left bdy\t"<<elementbdys[0]<<"\n";
	//cout <<"right bdy\t"<<elementbdys[nelts-1]<<"\n";
	if(x<elementbdys[0]){
		return 0;
	}
	if(x>elementbdys[nelts-1]){
		return nelts;
	}

	int lindx=0;
	int rindx=nelts-1;

	while((rindx-lindx)>1){
		int cindx=(lindx+rindx)/2;
		//cout << "lindx,cindx,rindx\t"<<lindx<<"\t"<<cindx<<"\t"<<rindx<<"\n";
		//cout << "x\t"<<x<<"\n";
		//cout << "lval,cval,rval\t"<<elementbdys[lindx]<<"\t"<<elementbdys[cindx]<<"\t"<<elementbdys[rindx]<<"\n";
		if(x<elementbdys[cindx]){
			rindx=cindx;
		}
		if(x>elementbdys[cindx]){
			lindx=cindx;
		}
	}
	//cout <<"returning lindx "<<lindx<<"\n";
	return lindx+1;
}
