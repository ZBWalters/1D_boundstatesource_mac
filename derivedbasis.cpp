#include "./classes.h"
//#include "./potential.h"
//#include "./basis.h"
//#include "./legendrebasis.h"
//#include "./laguerrebasis.h"
#include "./derivedbasis.h"

derivedbasis::derivedbasis(){
  order=0;
  parentorder=0;
  parentbasis=0;
  parentorder=0;
  omat=0;
  nablamat=0;
  nablasqmat=0;
  bfmat=0;
  Vmat=0;
  Hmat=0;
  xmin=0.;
  xmax=0.;
  ecsflag=false;
}

derivedbasis::derivedbasis(basis* bs_in, cmplx* bfmat_in, 
			   int neworder){
  //construct a new basis which is a linear combination of an existing basis
  //(The new basis may or may not have the same order as the underlying basis,
  //as it may prove convenient to throw out certain basis functions for the
  //purpose of boundary conditions, etc).

  //This is the generic constructor.  Specific constructors will be called if the parent is a legendre basis, laguerre basis, or another derived basis.

  parentbasis=bs_in;
  bfmat=arraycopy(bfmat_in,parentbasis->order*neworder);
  order=neworder;

  parentorder=parentbasis->order;

  //cout << "parent order "<<parentorder<<"\n";
  //cout << "order "<<order<<"\n";
  //cout <<"parentbasis->omat\n";
  //printmat(parentbasis->omat,parentorder,parentorder);
  //cout << "bfmat\n";
  //printmat(bfmat,order,neworder);
  omat=matrix_basischange(parentorder,order,parentbasis->omat,bfmat);
  cout << "derived omat\n";
  printomat();
  //printmat(omat,order,order);
  nablamat=matrix_basischange(parentorder,order,
			      parentbasis->nablamat,bfmat);
  //cout << "derived nablamat\n";
  //printmat(nablamat,order,order);
  nablasqmat=matrix_basischange(parentorder,order,
				parentbasis->nablasqmat,bfmat);
  //cout << "derived nablasqmat\n";
  //printmat(nablasqmat,order,order);

  //account for effect of changing scales on omat, nablasqmat (for nablamat,
  //scale factors cancel)
  cmplx scalefactor=abs(xmax-xmin)/2.;
  cmplx invscalefactor=1./scalefactor;
  arrayscale(omat,order*order,scalefactor);
  arrayscale(nablasqmat,order*order,invscalefactor);


  Vmat=initializearray(order*order);
  Hmat=initializearray(order*order);
  ecsflag=false;
}

derivedbasis::~derivedbasis(){
  parentbasis=0;
  //delete parentbasis;
  delete [] bfmat;
  delete [] omat;
  delete [] nablamat;
  delete [] nablasqmat;
  delete [] Vmat;
  delete [] Hmat;
  
	     
}

derivedbasis::derivedbasis(legendrebasis* bs_in, cmplx* bfmat_in, 
			   int neworder, rl xmin_in, rl xmax_in){
  //cout << "inside derived basis constructor\n";
  parentbasis=bs_in;
  bfmat=arraycopy(bfmat_in,parentbasis->order*neworder);
  order=neworder; 
  parentorder=parentbasis->order;
  omat=matrix_basischange(parentorder,order,parentbasis->omat,bfmat);
  nablamat=matrix_basischange(parentorder,order,
			      parentbasis->nablamat,bfmat);
  nablasqmat=matrix_basischange(parentorder,order,
				parentbasis->nablasqmat,bfmat);
  Vmat=initializearray(order*order);
  Hmat=initializearray(order*order);

  //cout << "omat intermediate \n";
  //printomat();
  xmin=xmin_in;
  xmax=xmax_in;

  //account for effect of changing scales on omat, nablasqmat (for nablamat,
  //scale factors cancel)
  cmplx scalefactor=(xmax-xmin)/2.;
  cmplx invscalefactor=1./scalefactor;
  arrayscale(omat,order*order,scalefactor);
  arrayscale(nablasqmat,order*order,invscalefactor);
  ecsflag=false;
  //cout << "omat after scaling \n";
  //printomat();

  //xtoy=[&](rl x) {map2intervals(x,xmin_in,xmax_in,-1.,1.);};
  //ytox=[&](rl y) {map2intervals(y,-1.,1.,xmin_in,xmax_in);};
  //cout << "order, parent order inside constructor "<<order<<" "<<parentorder<<"\n";
}

derivedbasis::derivedbasis(laguerrebasis* bs_in, cmplx* bfmat_in, 
			   int neworder, rl Rlag){
  parentbasis=bs_in;
  bfmat=arraycopy(bfmat_in,parentbasis->order*neworder);
  order=neworder; 
  parentorder=parentbasis->order;
  omat=matrix_basischange(parentorder,order,parentbasis->omat,bfmat);
  nablamat=matrix_basischange(parentorder,order,
			      parentbasis->nablamat,bfmat);
  nablasqmat=matrix_basischange(parentorder,order,
				parentbasis->nablasqmat,bfmat);
  Vmat=initializearray(order*order);
  Hmat=initializearray(order*order);
  //cout << "derived from laguerrebasis, kappa= "<<bs_in->kappa<<"\n";
  if(bs_in->kappa<0.){
    xmin=Rlag;
    xmax=0;
  }
  else{
    xmax=Rlag;
    xmin=0;
  }
  ecsflag=false;
}

cmplx* derivedbasis::evalbasisfuncs(rl x){
  //cout << "evalbasisfuncs order, parentorder "<<order<<" "<<parentorder<<"\n";
  cmplx* parentvals=parentbasis->evalbasisfuncs(parentbasis->mapfromchildrange(x,xmin,xmax));
  //cout << "parentvals "<<x<<", "<<parentbasis->mapfromchildrange(x,xmin,xmax)<<"\n";//xtoy(x,parentbasis)<<"\n";
  //printmat(parentvals,order,1);
  cmplx* retvals=vector_basischange(parentorder,order,parentvals,bfmat);
  delete [] parentvals;
  //cout << "basis function vals\n";
  //printmat(retvals,order,1);

  return retvals;
}

cmplx derivedbasis::evalcoeffvector(cmplx* coeffvector, rl x){
  //given a vector of basis function coeffs, evaluate the function
  //f(x)=\sum_{i} coeffvector_{i} X_{i}(x)
  //coeffvector is assumed to have length equal to order

  cmplx* bfvals=evalbasisfuncs(x);
  cmplx* retvalarray=initializearray(1);
  cmplx alpha=1.;
  cmplx beta=0.;
  cblas_zgemv(CblasColMajor,CblasNoTrans,1,order,&alpha,bfvals,1,
	      coeffvector,1,&beta,retvalarray,1);
  cmplx retval=retvalarray[0];
  delete [] retvalarray;
  delete [] bfvals;
  return retval;
}

str derivedbasis::printstr(cmplx* coeffvector){
  return printstr(order,coeffvector);
}

str derivedbasis::printstr(int printorder,cmplx* coeffvector){
  rl** inttable=integrationtable(printorder);
  rl* __restrict printpts=inttable[0];
  //printpts[0]=xmin;
  //printpts[printorder-1]=xmax;

  str retstr="";
  for(int i=0;i<printorder;i++){
    cmplx val=evalcoeffvector(coeffvector,printpts[i]);
    retstr+=to_string(printpts[i])+"\t"+to_string(real(val))+"\t"+
      to_string(imag(val))+"\n";
  }
  delete [] inttable[0];
  delete [] inttable[1];
  delete [] inttable;
  return retstr;
}

//void derivedbasis::Vmatsetup(potential* pot){
//  cmplx* prntVmat=parentbasis->Vmatsetup(potential* pot,xtoy,this);//parentbasis->funcmatrix(pot->V,xtoy,this);
//  Vmat=matrix_basischange(parentbasis->order,order,prntVmat,bfmat);
//  delete [] prntVmat;
//}

void derivedbasis::updateHmat(potential* pot){
  Vmatsetup(pot);
  if(Hmat!=0){
    delete [] Hmat;
  }
  Hmat=hmatsetup();
}

cmplx* derivedbasis::integrateVmat_testconvergence(potential* pot, rl x1, 
						   rl x2,  rl tol){
  cmplx* Vmat1=integrateVmat(pot,x1,x2);
  //cmplx* Vmat2=arraycopy(Vmat1,order*order);
  cmplx* Vmat2=integrateVmat(pot,x1,x2,2);
  //cout <<"integrateVmat_testconvergence "<<x1<<",\t"<<x2<<"\n";
  //cout<<"Vmat1\n";
  //printmat(Vmat1,order,order);
  //cout << "Vmat2\n";
  //printmat(Vmat2,order,order);
  cmplx* retmat=0;

  rl maxerr=0;
  for(int i=0;i<order*order;i++){
    cmplx avg=(Vmat1[i]+Vmat2[i])/2.;
    cmplx diff=Vmat2[i]-Vmat1[i];
    rl err=abs(diff);
    maxerr=max(maxerr,err);
  }
  if(maxerr<tol){
    delete [] Vmat1;
    retmat=Vmat2;
  }
  else{
    delete [] Vmat1;
    delete [] Vmat2;
    rl xavg=(x1+x2)/2.;
    cmplx* Vmatl=integrateVmat_testconvergence(pot,x1,xavg,tol);
    cmplx* Vmatr=integrateVmat_testconvergence(pot,xavg,x2,tol);
    retmat=new cmplx[order*order];
    for(int i=0;i<order*order;i++){
      retmat[i]=Vmatl[i]+Vmatr[i];
    }
    delete [] Vmatl;
    delete [] Vmatr;
  }
  return retmat;
}

cmplx* derivedbasis::integrateVmat(potential* pot, rl x1, 
					      rl x2, int nintervals){
//  cmplx* __restrict retmat=initializearray(order*order);
//  rl xavg=(x2+x1)/2.;
//  cmplx* __restrict V1mat=integrateVmat(pot,x1,xavg);
//  cmplx* __restrict V2mat=integrateVmat(pot,xavg,x2);
//  //cout << "Vmat left\n";
//  //printmat(V1mat,order,order);
//  //cout << "Vmat right\n";
//  //printmat(V2mat,order,order);
//  for(int i=0;i<order*order;i++){
//    retmat[i]=V1mat[i]+V2mat[i];
//  }
//  delete [] V1mat;
//  delete [] V2mat;
//  return retmat;
  cmplx* retmat=initializearray(order*order);
  rl dx=(x2-x1)/nintervals;
  for(int i=0;i<nintervals;i++){
    rl x1tmp=x1+dx*i;
    rl x2tmp=x1+dx*(i+1);
    cmplx* __restrict tmpVint=integrateVmat(pot,x1tmp,x2tmp);
    //cout << "tmpVint\n";
    //printmat(tmpVint,order,order);
    for(int j=0;j<order*order;j++){
      retmat[j]+=tmpVint[j];
    }
    delete [] tmpVint;
  }
  return retmat;
}


cmplx* derivedbasis::integrateVmat(potential* pot, rl x1, rl x2){
  //for a singular potential (or nearly so), it may prove desireable to have a
  //low order integrator rather than a high order one.  This routine uses the
  //midpoint method.
  rl xavg=(x1+x2)/2.;
  rl V=pot->V(xavg);
  rl dx=(x2-x1);
  cmplx* __restrict basisvals=this->evalbasisfuncs(xavg);
  cmplx* retmat=initializearray(order*order);
  for(int i=0;i<order;i++){
    for(int j=0;j<order;j++){
      retmat[arrayindx(i,j,order,order)]=basisvals[i]*basisvals[j]*V*dx;
    }
  }
  delete [] basisvals;
  return retmat;
}

cmplx* derivedbasis::integrateVmat(potential* pot, rl x1, rl x2, 
				   rl ecsbdy, rl ecstheta){
  //find matrix elements of Vecs, the analytic continuation of V
  cmplx* __restrict retmat=initializearray(order*order);
  int nint=2*order;
  rl** inttable=integrationtable(nint);
  rl* intpts=inttable[0];
  rl* intwts=inttable[1];
  for(int n=0;n<nint;n++){
    cmplx V=pot->Vecs(intpts[n],ecsbdy,ecstheta);
    cmplx* __restrict basisvals=this->evalbasisfuncs(intpts[n]);
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	retmat[arrayindx(i,j,order,order)]+=(basisvals[i]*basisvals[j]*
					     V*intwts[n]);
      }
    }
    delete [] basisvals;
  }
  delete [] intpts;
  delete [] intwts;
  delete [] inttable;

  //multiply entire array by exp(II*theta)
  cmplx II = cmplx(0.,1.);
  arrayscale(retmat,order*order,exp(-II*ecstheta));
  return retmat;
}

cmplx* derivedbasis::integratefuncmat(cmplx (*func)(rl x)){
  //find matrix elements of Vecs, the analytic continuation of V
  cmplx* __restrict retmat=initializearray(order*order);
  int nint=2*order;
  rl** inttable=integrationtable(nint);
  rl* intpts=inttable[0];
  rl* intwts=inttable[1];
  for(int n=0;n<nint;n++){
    cmplx Val=func(intpts[n]);
    cmplx* __restrict basisvals=this->evalbasisfuncs(intpts[n]);
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	retmat[arrayindx(i,j,order,order)]+=(basisvals[i]*basisvals[j]*
					     Val*intwts[n]);
      }
    }
    delete [] basisvals;
  }
  delete [] intpts;
  delete [] intwts;
  delete [] inttable;

  return retmat;
}

cmplx* derivedbasis::functimespsi(cmplx (*func)(rl x), cmplx* psi){
  cmplx* funcmat=integratefuncmat(func);
  cmplx* retvec=initializearray(order);
  cmplx alpha=1.;
  cmplx beta=0.;
  int inc=1;
  cblas_zgemv(CblasColMajor,CblasNoTrans,order,order,&alpha,
	      funcmat,order,psi,inc,&beta,retvec,inc);
  delete [] funcmat;
  return retvec;
}

//cmplx* derivedbasis::integrateVmat(potential* pot, rl x1, rl x2){
//  //perform integration using gauss legendre points and weights as this is
//  //quite high order, this may not be a desireable approach when the potential
//  //in question is singular or nearly so.  Then it might be better to use a
//  //lower order integration such as the midpoint method.
//  cmplx* retmat=initializearray(order*order);
//
//  int glorder=2*order;
//  rl* __restrict glpts=new rl[glorder];
//  rl* __restrict glwts=new rl[glorder];
//  gltable(glorder,glpts,glwts);
//  //map gauss legendre points & weights to current interval
//  for(int i=0;i<glorder;i++){
//    glwts[i]=glwts[i]*(x2-x1)/2.;
//    glpts[i]=x1+(x2-x1)*(glpts[i])/(2.);
//  }
//  cout << "glpts\n";
//  printmat(glpts,glorder,1);
//  cout << "glwts\n";
//  printmat(glwts,glorder,1);
//
//  for(int glk=0;glk<glorder;glk++){
//    rl V =pot->V(glpts[glk]);
//    cout << "V "<<glpts[glk]<<"\t"<<V<<"\n";
//    cmplx* __restrict basisvals=this->evalbasisfuncs(glpts[glk]);
//    cout << "basisvals\n";
//    printmat(basisvals,order,1);
//    for(int i=0;i<order;i++){
//      for(int j=0;j<order;j++){
//	//cout << "i,j,V"<<i<<" "<<j<<" "<<V<<"\n";
//	retmat[arrayindx(i,j,order,order)]+=basisvals[i]*basisvals[j]*
//	  V*glwts[glk];
//      }
//    }
//    delete [] basisvals;
//  }
//  delete [] glpts;
//  delete [] glwts;
//  cout << "retmat\n";
//  printmat(retmat,order,order);
//  return retmat;
//}

void derivedbasis::Vmatsetup(potential* pot){
  rl defaulttol=1.e-6;
  Vmatsetup(pot,defaulttol);
}

void derivedbasis::Vmatsetup(potential* pot, rl tol){
  if(Vmat!=0){
    delete [] Vmat;
  }
  Vmat=integrateVmat_testconvergence(pot,xmin,xmax,tol);
}

//void derivedbasis::Vmatsetup(potential* pot){
//  zeroarray(Vmat,order*order);
//  int glorder=2*order;
//  //cout << "calling integrationtable\n";
//  rl** inttable=integrationtable(glorder);
//  //cout << "back from integrationtable\n";
//  __restrict rl* intpts=inttable[0];
//  __restrict rl* intwts=inttable[1];
//  
//  //cout << "intpts\n";
//  //printmat(intpts,glorder,1);
//  //cout << "intwts\n";
//  //printmat(intwts,glorder,1);
//
//  
//  for(int glk=0;glk<glorder;glk++){
//    rl V =pot->V(intpts[glk]);
//    cmplx* __restrict basisvals=this->evalbasisfuncs(intpts[glk]);
//    //cout << "order"<<order<<"\n";
//    //cout << "basisvals\n";
//    //for(int i=0;i<order;i++){
//    //  cout <<"basisvals i "<<basisvals[i]<<"\n";
//    //}
//    for(int i=0;i<order;i++){
//      for(int j=0;j<order;j++){
//	//cout << "i,j,V"<<i<<" "<<j<<" "<<V<<"\n";
//	Vmat[arrayindx(i,j,order,order)]+=basisvals[i]*basisvals[j]*
//	  V*intwts[glk];
//	//cout << "i "<<i<<"\n";
//	//cout << "j"<<j<<"\n";
//	//cout <<"basisvals[i] "<<basisvals[i]<<"\n";
//	//cout <<"basisvals[j] "<<basisvals[j]<<"\n";
//	//cout <<"V "<<V<<"\n";
//	//cout <<"intwts[glk] "<<intwts[glk]<<"\n";
//	//cout << "Vmat[arrayindx(i,j,order,order)] "<<Vmat[arrayindx(i,j,order,order)]<<"\n";
//	//cout << "arrayindx(i,j,order,order) "<< arrayindx(i,j,order,order)<<"\n";
//      }
//    }
//    delete [] basisvals;
//  }
//
//  delete [] intpts;
//  delete [] intwts;
//  delete [] inttable;
//}

cmplx* derivedbasis::funcmatrix(cmplx (*func)(rl x), 
				rl (*ytox)(rl y, basis* bas1,
					   basis* bas2),basis* dbasis){
  return parentbasis->funcmatrix(func,ytox,dbasis);
}



cmplx* derivedbasis::hmatsetup(){
  //returns p^2/2m+Vmat
  cmplx* retmat=new cmplx[order*order];

  for(int i=0;i<order*order;i++){
    retmat[i]=-0.5*nablasqmat[i]+Vmat[i];
  }
  //cout << "hmat in hmatsetup\n";
  //printmat(retmat,order,order);
  return retmat;
}


//rl derivedbasis::xtoy_derived(rl x){
//  //xtoy if parent is a derived basis
//  return map2intervals(x,xmin,xmax,parentbasis->xmin,parentbasis->xmax);
//}
//
//rl derivedbasis::xtoy_legendre(rl x){
//  //xtoy if parent is a legendre basis (so that ymin=-1, ymax=1)
//  return map2intervals(x,xmin,xmax,-1.,1.);
//}
//
////rl derivedbasis::xtoy_laguerre(rl x){
////  return 
////}
//
//rl derivedbasis::ytox_derived(rl x){
//  //ytox if parent is a derived basis
//  return map2intervals(x,parentbasis->xmin,parentbasis->xmax,xmin,xmax);
//}
//
//rl derivedbasis::ytox_legendre(rl x){
//  //ytox if parent is a legendre basis (so that ymin=-1, ymax=1)
//  return map2intervals(x,-1.,1.,xmin,xmax);
//}
//
////rl derivedbasis::ytox_laguerre(rl x){
////  return 
////}

rl** derivedbasis::integrationtable(int glorder){
  //returns table of integration points, weights.  This is accomplished by
  //first calling integrationtable for the parent basis, then mapping the
  //points & weights to the current basis.  That way, the correct integration
  //table is returned regardless of whether the underlying basis consists of
  //legendre polynomials, laguerre polynomials, or some other basis.

  //cout << "calling parent integration table\n";
  //cout << "parent order "<<parentbasis->order<<"\n";
  rl** prntintegrationtable=parentbasis->integrationtable(glorder);
  //cout << "back from parent integration table\n";
  __restrict rl* prntintpts=prntintegrationtable[0];
  __restrict rl* prntintwts=prntintegrationtable[1];

  __restrict rl* intpts=new rl[glorder];
  __restrict rl* intwts=new rl[glorder];

  for(int i=0;i<glorder;i++){
    intpts[i]=parentbasis->mapy2range(prntintpts[i],xmin,xmax);//ytox(prntintpts[i],parentbasis);
  }

  rl dxdy=parentbasis->dy2dx(xmin,xmax);//dytodx(parentbasis);
  for(int i=0;i<glorder;i++){
    intwts[i]=prntintwts[i]*dxdy;
  }

  //cout << "xmin, xmax "<<xmin<<", "<<xmax<<"\n";
  //cout << "prntintpts\n";
  //printmat(prntintpts,order,1);
  //cout << "intpts\n";
  //printmat(intpts,order,1);

  delete [] prntintpts;
  delete [] prntintwts;
  delete [] prntintegrationtable;

  rl** rettable=new rl*[2];
  rettable[0]=intpts;
  rettable[1]=intwts;

  return rettable;
}

//methods related to change of basis
rl derivedbasis::mapy2range(rl y, rl x1, rl x2){
  rl xret=x1+(y-xmin)/(xmax-xmin)*(x2-x1);
  return xret;
}

rl derivedbasis::mapfromchildrange(rl childx, rl childxmin, rl childxmax){
  rl yret=xmin+(xmax-xmin)*(childx-childxmin)/(childxmax-childxmin);
  return yret;
}

rl derivedbasis::maptochildrange(rl y, rl childxmin, rl childxmax){
  rl xret=childxmin+(childxmax-childxmin)*(y-xmin)/(xmax-xmin);
  return xret;
}

rl derivedbasis::dy2dx(rl x1, rl x2){
  return (x2-x1)/(xmax-xmin);
}

//convert from legendre basis
rl derivedbasis::ytox(rl y, legendrebasis* lbas){
  // map parent basis to derived basis
  rl ymin=-1.;
  rl ymax=1.;
  
  rl xret=xmin+((y-ymin)/(ymax-ymin))*(xmax-xmin);
  return xret;
}

rl derivedbasis::xtoy(rl x, legendrebasis* lbas){
  rl ymin=-1.;
  rl ymax=1.;

  rl yret=ymin+((x-xmin)/(xmax-xmin))*(ymax-ymin);
  cout <<"xtoy function "<<x<<", "<<yret<<"\n";
  return yret;
}


rl derivedbasis::dytodx(legendrebasis* lbas){
  return (xmax-xmin)/2.;
}




//convert from derived basis
rl derivedbasis::ytox(rl y, derivedbasis* dbas1){
  rl ymin=dbas1->xmin;
  rl ymax=dbas1->xmax;

  rl xret=xmin+((y-ymin)/(ymax-ymin))*(xmax-xmin);
  return xret;
}

rl derivedbasis::xtoy(rl x, derivedbasis* dbas1){
  rl ymin=dbas1->xmin;
  rl ymax=dbas1->xmax;

  rl yret=ymin+((x-xmin)/(xmax-xmin))*(ymax-ymin);
  return yret;
}


rl derivedbasis::dytodx(derivedbasis* dbas1){
  rl ymin=dbas1->xmin;
  rl ymax=dbas1->xmax;
  return (xmax-xmin)/(ymax-ymin);
}


//convert from laguerre basis
rl derivedbasis::ytox(rl y, laguerrebasis* lagbas){
  if(lagbas->kappa<0){
    return xmin+y;
  }
  return xmax+y;
  
}

rl derivedbasis::xtoy(rl x, laguerrebasis* lagbas){
  if(lagbas->kappa<0){
    return x-xmin;
  }
  return x-xmax;
}


rl derivedbasis:: dytodx(laguerrebasis* lagbas){
  return 1.;
}

//generic basis change (y=x)
rl derivedbasis::ytox(rl y, basis* bas){
  return y;//default case is to do nothing
}

rl derivedbasis::xtoy(rl x, basis* bas){
  return x;//default case is to do nothing
}

rl derivedbasis::dytodx(basis* bas){
  return 1.;
}

rl derivedbasis::norm(cmplx* psi){
  return innerproduct(psi,psi);
}

rl derivedbasis::innerproduct(cmplx* psi1, cmplx* psi2){
  cmplx* __restrict psitmp=new cmplx[order];
  cmplx alpha=cmplx(1.,0.);
  cmplx beta=0.;
  cblas_zgemv(CblasColMajor,CblasNoTrans,order,order,&alpha,omat,order,psi2,1,
	      &beta,psitmp,1);
  //cmplx overlap=0.;
  //for(int i=0;i<order;i++){
  //  overlap+=conj(psi[i])*psitmp[i];
  //}
  //cout << "overlap sum "<<overlap<<"\n";

  //cout << "norm psi\n";
  //printmat(psi,order,1);
  //cout << "norm psitmp\n";
  //printmat(psitmp,order,1);
  
  cmplx* dot=new cmplx[1];
  cblas_zgemv(CblasColMajor,CblasConjTrans,order,1,&alpha,psi1,order,psitmp,1,
  	      &beta,dot,1);
  cmplx overlap=dot[0];
  //cout << "overlap zgemv "<<overlap<<"\n";
  
  delete [] psitmp;
  delete [] dot;
  return abs(overlap);
}

cmplx* derivedbasis::Olappsi(cmplx* psi_in){
  //returns Omat.psi
  return Xmatpsi(omat,psi_in);
}

cmplx* derivedbasis::Hpsi(cmplx* psi_in){
  //returns Omat.psi
  return Xmatpsi(Hmat,psi_in);
}

cmplx* derivedbasis::Vpsi(cmplx* psi_in){
  //returns Omat.psi
  return Xmatpsi(Vmat,psi_in);
}


cmplx* derivedbasis::Xmatpsi(cmplx* Xmat,cmplx* psi_in){
  //returns Omat.psi
  cmplx alpha=cmplx(1.,0.);
  cmplx beta=0.;
  cmplx* __restrict psiret=new cmplx[order];
  cblas_zgemv(CblasColMajor,CblasNoTrans,order,order,&alpha,Xmat,order,
	      psi_in,1,&beta,psiret,1);
  return psiret;
}

void derivedbasis::ecsrotate(rl theta){
    //rotate nabla matrix by exp(I theta), nablasq by exp(2 I theta)
  cmplx II =cmplx(0.,1.);

  //cout << "ecsrotate\t"<<theta<<"\n";

  for(int i=0;i<order*order;i++){
    ////omat[i]*=exp(II*theta);
    //nablamat[i]*=exp(II*theta);
    //nablasqmat[i]*=exp(2.*II*theta);
    

    //Subtle but important point: since omat and hmat will be added to omat
    //and hmat from other elements, you can't just multiply nablasqmat by
    //exp(2 II theta).  Instead, you have to multiply omat by exp(-II*theta)
    //and nablasqmat by exp(II*theta)=exp(-II*theta)*exp(2*II*theta)
    //(nablamat's factors of exp(-II*theta)*exp(II*theta) cancel.  This way,
    //the border functions which are nonzero in both the scaled and unscaled
    //regions are treated correctly.
    
    

    omat[i]*=exp(-II*theta);
    nablasqmat[i]*=exp(II*theta);
    Vmat[i]*=exp(-II*theta);
  }
  ecsflag=true;
  //for(int i=0;i<order*order;i++){
  //  Vmat[i]*=exp(-II*theta);//exp(-II*theta);//=0.;//zero potential matrix elements
  //}

}


void derivedbasis::ecsrotate(rl theta,potential* pot,rl ecsbdy,rl ecstheta){
  //rotate nabla matrix by exp(I theta), nablasq by exp(2 I theta)
  cmplx II =cmplx(0.,1.);

  //cout << "ecsrotate\t"<<theta<<"\n";

  for(int i=0;i<order*order;i++){
    ////omat[i]*=exp(II*theta);
    //nablamat[i]*=exp(II*theta);
    //nablasqmat[i]*=exp(2.*II*theta);
    

    //Subtle but important point: since omat and hmat will be added to omat
    //and hmat from other elements, you can't just multiply nablasqmat by
    //exp(2 II theta).  Instead, you have to multiply omat by exp(-II*theta)
    //and nablasqmat by exp(II*theta)=exp(-II*theta)*exp(2*II*theta)
    //(nablamat's factors of exp(-II*theta)*exp(II*theta) cancel.  This way,
    //the border functions which are nonzero in both the scaled and unscaled
    //regions are treated correctly.
    
    

    omat[i]*=exp(-II*theta);
    nablasqmat[i]*=exp(II*theta);
  }
  if(Vmat!=0){
    delete [] Vmat;
  }
  Vmat=integrateVmat(pot,xmin,xmax,ecsbdy,ecstheta);
  //for(int i=0;i<order*order;i++){
  //  Vmat[i]*=0.;//exp(-II*theta);//exp(-II*theta);//=0.;//zero potential matrix elements
  //}
  ecsflag=true;
}

cmplx* derivedbasis::funcolap(cmplx (*func)(rl x)){
//first, integrate function times basis set
  rl nintpts=2*order;
  rl** inttable=integrationtable(nintpts);
  rl* intpts=inttable[0];
  rl* intwts=inttable[1];
  //cout << "intpts\n";
  //printmat(intpts,order,1);
  //cout << "intwts\n";
  //printmat(intwts,order,1);
  cmplx* integrals=initializearray(order);
  for(int i=0;i<nintpts;i++){
    rl x=intpts[i];
    rl wt=intwts[i];
    cmplx funcval=func(x);
    cmplx* bfvals=evalbasisfuncs(x);
    for(int j=0;j<order;j++){
      integrals[j]+=funcval*bfvals[j]*wt;
    }
    delete [] bfvals;
  }

  delete [] intpts;
  delete [] intwts;
  delete [] inttable;
  return integrals;
}

cmplx* derivedbasis::funcolap(cmplx (*func)(rl x, rl y), rl param){
  //expand a function in terms of the basis set

  //first, integrate function times basis set
  rl** inttable=integrationtable(order);
  rl* intpts=inttable[0];
  rl* intwts=inttable[1];
  //cout << "intpts\n";
  //printmat(intpts,order,1);
  //cout << "intwts\n";
  //printmat(intwts,order,1);
  cmplx* integrals=initializearray(order);
  for(int i=0;i<order;i++){
    rl x=intpts[i];
    rl wt=intwts[i];
    cmplx funcval=func(x,param);
    cmplx* bfvals=evalbasisfuncs(x);
    for(int j=0;j<order;j++){
      integrals[j]+=funcval*bfvals[j]*wt;
    }
    delete [] bfvals;
  }

  delete [] intpts;
  delete [] intwts;
  delete [] inttable;
  return integrals;
}



cmplx* derivedbasis::funcvector(cmplx (*func)(rl x)){
  //expand a function in terms of the basis set

  cmplx* integrals=funcolap(func);

  //cout << "integrals before solve\n";
  //printmat(integrals,order,1);
  //cout << "omat before solve\n";
  //printmat(omat,order,order);
  //cout << "parent omat\n";
  //printmat(parentbasis->omat,order,order);

  //now solve linear system for basis function coefficients
  cmplx* lhsmat=arraycopy(omat,order*order);
  int nrhs=1;
  int* ipiv=new int[order];
  int info=0;
  //zgesv_(&order,&nrhs,lhsmat,&order,ipiv,integrals,&order,&info);
  //if(info!=0){
  //  cout << "info ne 0 after zgesv_ in funcvector!";
  //}
  char uplo='U';
  cmplx* tmpwork=new cmplx[1];
  int lwork=-1;
  zhesv_(&uplo, &order, &nrhs, reinterpret_cast <Lcmplx*> (lhsmat), &order,
		  ipiv, reinterpret_cast <Lcmplx*> (integrals), &order,
		  reinterpret_cast <Lcmplx*> (tmpwork), &lwork, &info);
  if(info!=0){
    cout << "problem with zhesv_ in funcvector 1st call!\n";
  }
  lwork=int(real(tmpwork[0]));
  delete [] tmpwork;
  cmplx* work=new cmplx[lwork];
  zhesv_(&uplo, &order, &nrhs, reinterpret_cast <Lcmplx*> (lhsmat), &order,
		  ipiv, reinterpret_cast <Lcmplx*> (integrals), &order,
		  reinterpret_cast <Lcmplx*> (work), &lwork, &info);
  if(info!=0){
    cout << "problem with zhesv_ in funcvector 1st call!\n";
  }
  
  delete [] ipiv;
  delete [] lhsmat;

  return integrals;
}


cmplx* derivedbasis::funcvector(cmplx (*func)(rl x, rl y), rl param){
  //expand a function in terms of the basis set

  //cout << "intpts\n";
  //printmat(intpts,order,1);
  //cout << "intwts\n";
  //printmat(intwts,order,1);
  //first, integrate function against basis set
  cmplx* integrals=funcolap(func,param);
  
  //cout << "integrals before solve\n";
  //printmat(integrals,order,1);
  //cout << "omat before solve\n";
  //printmat(omat,order,order);
  //cout << "parent omat\n";
  //printmat(parentbasis->omat,order,order);

  //now solve linear system for basis function coefficients
  cmplx* lhsmat=arraycopy(omat,order*order);
  int nrhs=1;
  int* ipiv=new int[order];
  int info=0;
  zgesv_(&order, &nrhs, reinterpret_cast <Lcmplx*> (lhsmat), &order, ipiv,
		  reinterpret_cast <Lcmplx*> (integrals), &order, &info);
  if(info!=0){
    cout << "info ne 0 after zgesv_ in funcvector!";
  }
  //cout << "integrals\n";
  //printmat(integrals,order,1);
  delete [] ipiv;
  delete [] lhsmat;
  return integrals;
}

