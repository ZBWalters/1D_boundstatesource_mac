#include "./classes.h"
#include "./basis.h"
#include "./laguerrebasis.h"

//The tricky detail is that the gauss laguerre weights
//already include the factor of exp(-x) needed to give orthogonality, so that
//\sum_{i} L_n(x_i) L_m(x_i) w_i = delta_{n,m} In the code, I would like to use
//basis functions consisting of Laguerre polynomials times exponentials, or
//bf_n(x)=L_n(x) exp(-x/2), so that \int_{0}^{\infty} bf_n(x)
//bf_m(x)=delta_{n,m}.  With respect to this basis, expansion coefficients are
//found via
// A_n=\sum_i f(x_i) exp(x_i/2) bf_n(x_i) w_i



laguerrebasis::laguerrebasis(int order_in, rl kappa_in){
  order=order_in;
  kappa=kappa_in;
  omat=initializearray(order*order);
  nablamat=initializearray(order*order);
  nablasqmat=initializearray(order*order);


  int glorder=2*order;
  rl* glpts=new rl[glorder];
  rl* glwts=new rl[glorder];

  //overlap matrix
  laguerretable(glorder,glpts,glwts);
  for(int glk=0;glk<glorder;glk++){
    rl* lvals=laguerrevalues(glorder,glpts[glk]);
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	omat[arrayindx(i,j,order,order)]+=lvals[i]*lvals[j]*glwts[glk];
      }
    }
    delete [] lvals;
  }
  //cout << "laguerre omat\n";
  //printmat(omat,order,order);

  //nabla matrix (first derivs)
  for(int glk=0;glk<glorder;glk++){
    rl* lvals=laguerrevalues(glorder,glpts[glk]);
    rl* lpvals=laguerrederivs(glorder,0,1,glpts[glk]);
    //cout << "laguerre point "<<glpts[glk]<<"\n";
    //cout << "laguerre values\n";
    //printmat(lvals,glorder,1);
    //cout << "laguerre derivs\n";
    //printmat(lpvals,glorder,1);
    
 //   if(kappa >0){//kappa >0 implies that x is less than zero
 //     for(int i=0;i<glorder;i++){
 //	lpvals[i]*=-1.;
 //     }
 //   }

    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	nablamat[arrayindx(i,j,order,order)]
	  +=lvals[i]*(lpvals[j]-0.5*lvals[j])*glwts[glk];
      }
    }
    delete [] lvals;
    delete [] lpvals;
  }
  //if(kappa>0){//kappa >0 implies that x is less than zero
  //arrayscale(nablamat,order*order,0.);//test if p term is hurting
  //}

  //cout << "laguerre nablamat\n";
  //printmat(nablamat,order,order);



  //nablasq matrix (2nd derivs)
   for(int glk=0;glk<glorder;glk++){
    rl* lvals=laguerrevalues(glorder,glpts[glk]);
    rl* lpvals=laguerrederivs(glorder,0,1,glpts[glk]);
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	nablasqmat[arrayindx(i,j,order,order)]-=(lpvals[i]-0.5*lvals[i])*
	  (lpvals[j]-0.5*lvals[j])*glwts[glk];
	
      }
    }
    delete [] lvals;
    delete [] lpvals;
   }
   

   //finally, have to adjust nabla and nablasq matrices to account for nonzero
   //values of kappa
   //recall that kappa >0 implies that x is less than zero
   arrayscale(nablamat,order*order,-1./kappa);
   arrayscale(nablasqmat,order*order,pow(kappa,-2.));

   cout << "kappa\t"<<kappa<<"\n";
   cout << "laguerre nablasqmat\n";
   printmat(nablasqmat,order,order);

   //set nablamat to zero so that vector potential doesn't cause divergence
   //arrayscale(nablamat,order*order,0.);


  
  
//  //nabla matrix (first derivs)
//  for(int glk=0;glk<glorder;glk++){
//    rl* lvals=laguerrevalues(glorder,glpts[glk]);
//    rl* lpvals=laguerrederivs(glorder,0,1,glpts[glk]);
//    if(kappa >0){//kappa >0 implies that x is less than zero
//      for(int i=0;i<glorder;i++){
//	lpvals[i]*=-1.;
//      }
//    }
//    for(int i=0;i<order;i++){
//      for(int j=0;j<order;j++){
//	nablamat[arrayindx(i,j,order,order)]+=lvals[i]*lpvals[j]*
//	  exp(-abs(kappa*glpts[glk]))*glwts[glk];
//      }
//    }
//    delete [] lvals;
//    delete [] lpvals;
//  }
//
//  //nablasq matrix (2nd derivs)
//  for(int glk=0;glk<glorder;glk++){
//    rl* lvals=laguerrevalues(glorder,glpts[glk]);
//    rl* lpvals=laguerrederivs(glorder,0,2,glpts[glk]);
//    for(int i=0;i<order;i++){
//      for(int j=0;j<order;j++){
//	nablasqmat[arrayindx(i,j,order,order)]+=lvals[i]*lpvals[j]*
//	  exp(-abs(kappa*glpts[glk]))*glwts[glk];
//      }
//    }
//    delete [] lvals;
//    delete [] lpvals;
//  }
  
  delete [] glpts;
  delete [] glwts;

}

laguerrebasis::~laguerrebasis(){
  delete [] omat;
  delete [] nablamat;
  delete [] nablasqmat;
}

cmplx* laguerrebasis::funcmatrix(cmplx (*func)(rl x), rl (*ytox)(rl y, basis* bas1, basis* bas2), basis* dbasis){
  cmplx* retmat=initializearray(order*order);


  int glorder=2*order;
  rl* glpts=new rl[glorder];
  rl* glwts=new rl[glorder];
  laguerretable(glorder,glpts,glwts);
  for(int glk=0;glk<glorder;glk++){
    rl* lvals=laguerrevalues(order,glpts[glk]);
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	retmat[arrayindx(i,j,order,order)]+=lvals[i]*lvals[j]*glwts[glk]*
	  func(ytox(glpts[glk],this,dbasis));
      }
    }
  }
  return retmat;
}

cmplx* laguerrebasis::evalbasisfuncs(rl y){
  //return L_n(2 kappa r)* exp(-kappa r), as in 
  //Scrinzi, PHYSICAL REVIEW A 81, 053845 (2010)

  rl tmpy=abs(kappa*y);
  cmplx* retvals=evallaguerrefuncs(tmpy);
  for(int i=0;i<order;i++){
    retvals[i]*=exp(-0.5*abs(kappa*y));
  }
  //cout << "laguerre values "<<y<<"\n";
  //printmat(retvals,order,1);

  return retvals;

}

cmplx* laguerrebasis::evallaguerrefuncs(rl y){
  rl* lvals=laguerrevalues(order,abs(y));
  cmplx* retvals=new cmplx[order];
  for(int i=0;i<order;i++){
    retvals[i]=lvals[i];
  }
  delete [] lvals;
  return retvals;
}

rl** laguerrebasis::integrationtable(int glorder){
  rl* glpts=new rl[glorder];
  rl* glwts=new rl[glorder];
  laguerretable(glorder,glpts,glwts);
  rl** rettable=new rl*[2];
  rettable[0]=glpts;
  rettable[1]=glwts;
  if(kappa>0){
    rl* negglpts=new rl[glorder];
    rl* negglwts=new rl[glorder];
    for(int i=0;i<glorder;i++){
      negglpts[i]=glpts[glorder-1-i];
      negglwts[i]=glwts[glorder-1-i]*exp(0.5*abs(kappa*negglpts[i]));
    }
    rettable[0]=negglpts;
    rettable[1]=negglwts;
    delete [] glpts;
    delete [] glwts;
  }
  return rettable;
}

rl laguerrebasis::mapy2range(rl y, rl x1, rl x2){
  //cout << "laguerre mapy2range\n";
  //cout << "kappa, y, x1, x2\n";
  //cout <<kappa<<" "<<y<<" "<<x1<<" "<<x2<<"\n";
  rl xret=0.;
  if(kappa>0.){
    xret=x2-y;
  }
  if(kappa<=0.){
    xret=x1+y;
  }
  //cout << "xret "<<xret<<"\n\n";
  return xret;
}

rl laguerrebasis::mapfromchildrange(rl childx, rl childxmin, rl childxmax){
  rl yret=0.;
  if(kappa>0.){
    yret=childxmax-childx;
  }
  if(kappa<=0.){
    yret=childx-childxmin;
  }
  return yret;
}

rl laguerrebasis::maptochildrange(rl y, rl childxmin, rl childxmax){
  rl xret=0.;
  if(kappa>0.){
    xret=childxmin+y;
  }
  if(kappa<=0.){
    xret=childxmax-y;
  }
  return xret;
}

rl laguerrebasis::dy2dx(rl x1, rl x2){
  return 1.;
}

//cmplx* laguerrebasis::xpsi(cmplx* __restrict psi){
//  //return 0, since psi will always be zero in this region for problems of
//  //interest
//  cmplx* __restrict retarray=initializerray(order);
//  return retarray;
//}
