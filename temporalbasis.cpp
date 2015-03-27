#include "./classes.h"
//#include "./basis.h"
//#include "./legendrebasis.h"
#include "./temporalbasis.h"

//temporalbasis::temporalbasis(derivedbasis* prnt_in,rl dt){
//  prnt=prnt_in;
//  order=prnt->order;
//  umat=new cmplx[order*order];
//  qmat=new cmplx[order*order];
//  Amat=initializearray(order*order);
//  bfmat=prnt->bfmat;
//  for(int i=0;i<order*order;i++){
//    umat[i]=prnt->omat[i]*dt;
//    qmat[i]=prnt->nablamat[i];
//  }
//}

temporalbasis::temporalbasis(legendrebasis* bas, cmplx* bfmat_in, 
			     rl tstart,rl tstop): tmin(xmin), tmax(xmax),umat(omat),qmat(nablamat) {
  //derivedbasis(bas,bfmat_in,bas->order,tstart,tstop);//have to set these terms directly
  parentbasis=bas;
  parentorder=bas->order;
  order=bas->order;
  bfmat=arraycopy(bfmat_in,order*parentorder);
  omat=matrix_basischange(parentorder,order,parentbasis->omat,bfmat);
  nablamat=matrix_basischange(parentorder,order,
			      parentbasis->nablamat,bfmat);

  ////test: should nablamat be transposed?
  //cout << "nablamat before transpose\n";
  //printmat(nablamat,order,order);
  //cmplx* tmpmat=nablamat;
  //nablamat=transpose(tmpmat,order,order);
  //delete [] tmpmat;
  //cout << "nablamat\n";
  //printmat(nablamat,order,order);
  //cout << "qmat\n";
  //printmat(qmat,order,order);

  nablasqmat=matrix_basischange(parentorder,order,
				parentbasis->nablasqmat,bfmat);
  Vmat=initializearray(order*order);

  //account for different range of temporal basis and underlying legendre basis 
  cmplx scalefactor=(tstop-tstart)/2.;
  cmplx invscalefactor=1./scalefactor;
  arrayscale(omat,order*order,scalefactor);
  arrayscale(nablasqmat,order*order,invscalefactor);

  //cout << "back from derived basis constructor "<<order<<" "<<parentorder<<"\n";
  //umat=omat;//this should be a reference, not just assignment
  //qmat=nablamat;//this should be a reference, not just assignment
  //rl& tmin=xmin;
  //rl& tmax=xmax;
  tmin=tstart;
  tmax=tstop;

  //cout << "tmin, xmin in constructor "<<tmin<<"\t"<<xmin<<"\n";
  //cout << "tmax, xmax in constructor "<<tmax<<"\t"<<xmax<<"\n";
  order=bas->order;
  parentorder=bas->order;
  //cout << "parent basis order "<<bas->order<<" "<<parentorder<<"\n";
  //cout << "temporal basis order "<<order<<"\n";
  Amat=initializearray(order*order);
  Asqmat=initializearray(order*order);
  Emat=initializearray(order*order);
  
}

temporalbasis::temporalbasis(legendrebasis* bas, cmplx* bfmat_in, pulse* Pls, 
			     rl tstart,rl tstop): tmin(xmin), tmax(xmax),umat(omat),qmat(nablamat) {
  //derivedbasis(bas,bfmat_in,bas->order,tstart,tstop);//have to set these terms directly
  parentbasis=bas;
  parentorder=bas->order;
  order=bas->order;
  bfmat=arraycopy(bfmat_in,order*parentorder);
  omat=matrix_basischange(parentorder,order,parentbasis->omat,bfmat);
  nablamat=matrix_basischange(parentorder,order,
			      parentbasis->nablamat,bfmat);

  ////test: should nablamat be transposed?
  //cout << "nablamat before transpose\n";
  //printmat(nablamat,order,order);
  //cmplx* tmpmat=nablamat;
  //nablamat=transpose(tmpmat,order,order);
  //delete [] tmpmat;
  //cout << "nablamat\n";
  //printmat(nablamat,order,order);
  //cout << "qmat\n";
  //printmat(qmat,order,order);

  nablasqmat=matrix_basischange(parentorder,order,
				parentbasis->nablasqmat,bfmat);
  Vmat=initializearray(order*order);

  //account for different range of temporal basis and underlying legendre basis 
  cmplx scalefactor=(tstop-tstart)/2.;
  cmplx invscalefactor=1./scalefactor;
  arrayscale(omat,order*order,scalefactor);
  arrayscale(nablasqmat,order*order,invscalefactor);

  //cout << "back from derived basis constructor "<<order<<" "<<parentorder<<"\n";
  //umat=omat;//this should be a reference, not just assignment
  //qmat=nablamat;//this should be a reference, not just assignment
  //rl& tmin=xmin;
  //rl& tmax=xmax;
  tmin=tstart;
  tmax=tstop;

  //cout << "tmin, xmin in constructor "<<tmin<<"\t"<<xmin<<"\n";
  //cout << "tmax, xmax in constructor "<<tmax<<"\t"<<xmax<<"\n";
  order=bas->order;
  parentorder=bas->order;
  //cout << "parent basis order "<<bas->order<<" "<<parentorder<<"\n";
  //cout << "temporal basis order "<<order<<"\n";
  Amat=initializearray(order*order);
  Asqmat=initializearray(order*order);
  Emat=initializearray(order*order);

  updateTmats(tstart,tstop,Pls);
}




temporalbasis::temporalbasis(temporalbasis* tbas, cmplx* bfmat_in, pulse* Pls, 
			     rl tstart,rl tstop): tmin(xmin), tmax(xmax),umat(omat),qmat(nablamat) {
  //derivedbasis(bas,bfmat_in,bas->order,tstart,tstop);//have to set these terms directly
  basis* bas=tbas->parentbasis;
  parentbasis=bas;
  parentorder=bas->order;
  order=bas->order;
  bfmat=arraycopy(bfmat_in,order*parentorder);
  omat=matrix_basischange(parentorder,order,parentbasis->omat,bfmat);
  nablamat=matrix_basischange(parentorder,order,
			      parentbasis->nablamat,bfmat);

  ////test: should nablamat be transposed?
  //cout << "nablamat before transpose\n";
  //printmat(nablamat,order,order);
  //cmplx* tmpmat=nablamat;
  //nablamat=transpose(tmpmat,order,order);
  //delete [] tmpmat;
  //cout << "nablamat\n";
  //printmat(nablamat,order,order);
  //cout << "qmat\n";
  //printmat(qmat,order,order);

  nablasqmat=matrix_basischange(parentorder,order,
				parentbasis->nablasqmat,bfmat);
  Vmat=initializearray(order*order);

  //account for different range of temporal basis and underlying legendre basis 
  cmplx scalefactor=(tstop-tstart)/2.;
  cmplx invscalefactor=1./scalefactor;
  arrayscale(omat,order*order,scalefactor);
  arrayscale(nablasqmat,order*order,invscalefactor);

  //cout << "back from derived basis constructor "<<order<<" "<<parentorder<<"\n";
  //umat=omat;//this should be a reference, not just assignment
  //qmat=nablamat;//this should be a reference, not just assignment
  //rl& tmin=xmin;
  //rl& tmax=xmax;
  tmin=tstart;
  tmax=tstop;

  //cout << "tmin, xmin in constructor "<<tmin<<"\t"<<xmin<<"\n";
  //cout << "tmax, xmax in constructor "<<tmax<<"\t"<<xmax<<"\n";
  order=bas->order;
  parentorder=bas->order;
  //cout << "parent basis order "<<bas->order<<" "<<parentorder<<"\n";
  //cout << "temporal basis order "<<order<<"\n";
  Amat=initializearray(order*order);
  Asqmat=initializearray(order*order);
  Emat=initializearray(order*order);

  updateTmats(tstart,tstop,Pls);
}











temporalbasis::~temporalbasis(){
  //don't delete prnt, because same basis might be used as the parent for
  //multiple temporal bases
  //delete [] umat;
  //delete [] qmat;
  delete [] Amat;
  delete [] Asqmat;
  delete [] Emat;
}

void temporalbasis::updateTmats(rl t1, rl t2, pulse* Pls){
  delete [] omat;
  delete [] nablamat;
  delete [] nablasqmat;
  tmin=t1;
  tmax=t2;
  omat=matrix_basischange(parentorder,order,parentbasis->omat,bfmat);
  nablamat=matrix_basischange(parentorder,order,
			      parentbasis->nablamat,bfmat);
  nablasqmat=matrix_basischange(parentorder,order,
				parentbasis->nablasqmat,bfmat);

  //account for different range of temporal basis and underlying legendre basis 
  cmplx scalefactor=(tmax-tmin)/2.;
  cmplx invscalefactor=1./scalefactor;
  arrayscale(omat,order*order,scalefactor);
  arrayscale(nablasqmat,order*order,invscalefactor);

  //update Amat
  Amatsetup(Pls);
}

void temporalbasis::Amatsetup(pulse* Pls){
  zeroarray(Amat,order*order);
  zeroarray(Asqmat,order*order);
  zeroarray(Emat,order*order);
  int glorder=2*order;
  //cout << "calling integrationtable\n";
  rl** inttable=integrationtable(glorder);
  //cout << "back from integrationtable\n";
  __restrict rl* intpts=inttable[0];
  __restrict rl* intwts=inttable[1];
  
  //cout << "intpts\n";
  //printmat(intpts,glorder,1);
  //cout << "intwts\n";
  //printmat(intwts,glorder,1);

  
  for(int glk=0;glk<glorder;glk++){
    rl A =Pls->Az(intpts[glk]);
    //cout <<"Amatsetup Az "<<A<<"\n";
    cmplx* __restrict basisvals=this->evalbasisfuncs(intpts[glk]);
    //cout << "order"<<order<<"\n";
    //cout << "basisvals\n";
    //for(int i=0;i<order;i++){
    //  cout <<"basisvals i "<<basisvals[i]<<"\n";
    //}
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	//cout << "i,j,V"<<i<<" "<<j<<" "<<V<<"\n";
	Amat[arrayindx(i,j,order,order)]+=basisvals[i]*basisvals[j]*
	  A*intwts[glk];
	//cout << "i "<<i<<"\n";
	//cout << "j"<<j<<"\n";
	//cout <<"basisvals[i] "<<basisvals[i]<<"\n";
	//cout <<"basisvals[j] "<<basisvals[j]<<"\n";
	//cout <<"V "<<V<<"\n";
	//cout <<"intwts[glk] "<<intwts[glk]<<"\n";
	//cout << "Vmat[arrayindx(i,j,order,order)] "<<Vmat[arrayindx(i,j,order,order)]<<"\n";
	//cout << "arrayindx(i,j,order,order) "<< arrayindx(i,j,order,order)<<"\n";
      }
    }
    delete [] basisvals;
  }

    for(int glk=0;glk<glorder;glk++){
    rl E =Pls->Ez(intpts[glk]);
    //cout <<"Amatsetup Az "<<A<<"\n";
    cmplx* __restrict basisvals=this->evalbasisfuncs(intpts[glk]);
    //cout << "order"<<order<<"\n";
    //cout << "basisvals\n";
    //for(int i=0;i<order;i++){
    //  cout <<"basisvals i "<<basisvals[i]<<"\n";
    //}
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	//cout << "i,j,V"<<i<<" "<<j<<" "<<V<<"\n";
	Emat[arrayindx(i,j,order,order)]+=basisvals[i]*basisvals[j]*
	  E*intwts[glk];
	//cout << "i "<<i<<"\n";
	//cout << "j"<<j<<"\n";
	//cout <<"basisvals[i] "<<basisvals[i]<<"\n";
	//cout <<"basisvals[j] "<<basisvals[j]<<"\n";
	//cout <<"V "<<V<<"\n";
	//cout <<"intwts[glk] "<<intwts[glk]<<"\n";
	//cout << "Vmat[arrayindx(i,j,order,order)] "<<Vmat[arrayindx(i,j,order,order)]<<"\n";
	//cout << "arrayindx(i,j,order,order) "<< arrayindx(i,j,order,order)<<"\n";
      }
    }
    delete [] basisvals;
  }
    //cout << "E mat\n";
    //printmat(Emat,order,order);

//    //new way to get Asqmat: square Amat
//    //first, find inverse of overlap matrix
//    cmplx* tmpu=arraycopy(umat,order*order);//copy overlap matrix
//    int* ipiv=new int[order];
//    int lwork=-1;
//    cmplx* tmpwork=new cmplx[1];
//    int info=0;
//    zgetrf(&order,&order,tmpu,&order,ipiv,&info);
//    if(info!=0){
//      cout << "problem with zgetrf!\n";
//    }
//    zgetri(&order,tmpu,&order,ipiv,tmpwork,&lwork,&info);
//    if(info!=0){
//      cout << "problem with zgetri 1!\n";
//    }
//    lwork=int(real(tmpwork[0]));
//    delete [] tmpwork;
//    //cout << "lwork "<<lwork<<"\n";
//    cmplx* work=new cmplx[lwork];
//    zgetri(&order,tmpu,&order,ipiv,work,&lwork,&info);
//    if(info!=0){
//      cout << "problem with zgetri 2!\n";
//    }
//    delete [] work;
//    delete [] ipiv;
//    //cout << "umat inverse\n";
//    //printmat(tmpu,order,order);
//    
//    cmplx alpha=1.;
//    cmplx beta=0.;
//    cmplx* UinvA=new cmplx[order*order];
//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
//		order,order,order,&alpha,tmpu,order,
//		Amat,order,&beta,UinvA,order);
//    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
//		order,order,order,&alpha,Amat,order,
//		UinvA,order,&beta,Asqmat,order);
//    delete [] UinvA;

    //old way to get Asqmat: integrate Asq wrt time
  for(int glk=0;glk<glorder;glk++){
    rl Asq =pow(Pls->Az(intpts[glk]),2.);
    //cout <<"Amatsetup Azsq "<<Asq<<"\n";
    cmplx* __restrict basisvals=this->evalbasisfuncs(intpts[glk]);
    //cout << "order"<<order<<"\n";
    //cout << "basisvals\n";
    //for(int i=0;i<order;i++){
    //  cout <<"basisvals i "<<basisvals[i]<<"\n";
    //}
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	//cout << "i,j,V"<<i<<" "<<j<<" "<<V<<"\n";
	Asqmat[arrayindx(i,j,order,order)]-=basisvals[i]*basisvals[j]*
	  Asq*intwts[glk];
	
	//cout << "i "<<i<<"\n";
	//cout << "j"<<j<<"\n";
	//cout <<"basisvals[i] "<<basisvals[i]<<"\n";
	//cout <<"basisvals[j] "<<basisvals[j]<<"\n";
	//cout <<"V "<<V<<"\n";
	//cout <<"intwts[glk] "<<intwts[glk]<<"\n";
	//cout << "Vmat[arrayindx(i,j,order,order)] "<<Vmat[arrayindx(i,j,order,order)]<<"\n";
	//cout << "arrayindx(i,j,order,order) "<< arrayindx(i,j,order,order)<<"\n";
      }
    }
    delete [] basisvals;
  }
 

  delete [] intpts;
  delete [] intwts;
  delete [] inttable;

}

cmplx* temporalbasis::Eiwtvec(rl En){
  //find expansion of exp[II*En*t] wrt temporal basis
  //ie, solve II O Q - H U=0
  cmplx II=cmplx(0.,1.);
  
  int redorder=order-1;

  //rhsvec is action residual resulting from initial condition 
  //f(tstart)=exp(-II*tstart)
  cmplx startval=exp(-II*En*tmin);
  cmplx* __restrict rhsvec=initializearray(redorder);
  for(int n=1;n<order;n++){
    rhsvec[n-1]=-startval*(II*qmat[arrayindx(n,0,order,order)]-
		 En*umat[arrayindx(n,0,order,order)]);
  }

  cmplx* __restrict lhsmat=initializearray((redorder)*(redorder));
  for(int n=1;n<order;n++){
    for(int m=1;m<order;m++){
      lhsmat[arrayindx(n-1,m-1,redorder,redorder)]=
	II*qmat[arrayindx(n,m,order,order)]-En*umat[arrayindx(n,m,order,order)];
    }
  }
    
  int* ipiv=new int[redorder];
  int info=0;
  int nrhs=1;
  zgesv_(&redorder, &nrhs, reinterpret_cast <Lcmplx*> (lhsmat), &redorder,
		  ipiv, reinterpret_cast <Lcmplx*> (rhsvec), &redorder, &info);
  if(info!=0){
    cout << "problem with zgesv_ in Expiwt!\n";
  }


  cmplx* __restrict retvec=initializearray(order);
  retvec[0]=startval;
  for(int n=1;n<order;n++){
    retvec[n]=rhsvec[n-1];
  }


    delete [] ipiv;
    delete [] lhsmat;
    delete [] rhsvec;
  
  return retvec;
}
