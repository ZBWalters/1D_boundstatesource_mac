#include "./classes.h"
#include "./basis.h"
#include "./legendrebasis.h"

//class legendrebasis : public basis;

//legendrebasis::legendrebasis(legendrebasis* inpbasis){
//  //copy omat,nablamat, nablasqmat from inpbasis
//  order=inpbasis->order;
//  omat=arraycopy(inpbasis->omat,order*order);
//  nablamat=arraycopy(inpbasis->nablamat,order*order);
//  nablasqmat=arraycopy(inpbasis->nablasqmat,order*order);
//}

legendrebasis::legendrebasis(int order_in){
  
  order=order_in;
  
  omat=initializearray(order*order);
  nablamat=initializearray(order*order);
  nablasqmat=initializearray(order*order);

  for (int i=0;i<order;i++){
    omat[arrayindx(i,i,order,order)]=2./(2.*i+1.);
  }

  for(int i=0;i<order;i++){
    for(int j=0;j<order;j++){
      //if((i>j) and ((i-j)%2==0)){
      //	nablasqmat[arrayindx(i,j,order,order)]=(i*(i+1.)-(pow(j,2.)+j));
      //	//psqmat[arrayindx(i,j,order,order)]/=2.;//KE=p^2/2
      //}
      //if((j>=i) and ((j-i)%2==0)){
      //	nablasqmat[arrayindx(i,j,order,order)]=(j*(j+1.)-(pow(i,2.)+i));
      //	//psqmat[arrayindx(i,j,order,order)]/=2.;//KE=p^2/2
      //}
      if((j-i)%2==0){
      int k=min(i,j);
      nablasqmat[arrayindx(i,j,order,order)]=-rl(k)*(rl(k)+1.);
      }
    }
  }

  
  for(int i=0;i<order;i++){
    for(int j=0;j<order;j++){
      //nablamat[arrayindx(i,j,order,order)]=0.;//since A is constant, its
      //spatial derivative is zero
      if((j>i) and ((j-i)%2==1)){
  	nablamat[arrayindx(i,j,order,order)]=2.;//2.;//2.;
  	//Factor of 2 comes from Legendre
  	//Polynomials
      }
    }
  }



//  cout<<"omat\n";
//  printmat(omat,order,order);
//
//  cout<<"nablamat\n";
//  printmat(nablamat,order,order);
//
//  cout<<"nablasqmat\n";
//  printmat(nablasqmat,order,order);
  
    
  
//  basis(order,omat,nablamat,nablasqmat);
//  delete [] omat;
//  delete [] nablamat;
//  delete [] nablasqmat;
}

legendrebasis::~legendrebasis(){
  delete [] omat;
  delete [] nablamat;
  delete [] nablasqmat;
}

cmplx* legendrebasis::funcmatrix(cmplx (*func)(rl x), rl (*ytox)(rl y, basis* bas1, basis* bas2),basis* dbasis){
  //return matrix of int_{-1}^{1}dy P_i(y)*func(x(y))*P_j(y) *dx/dy
  //ytox converts from units useful in legendre polynomials to units useful
  //in evaluating func.
  cmplx* retmat=new cmplx[order*order];
  zeroarray(retmat,order*order);

  rl dxdy=(ytox(1.,this,dbasis)-ytox(-1.,this,dbasis))/2.;

  int glorder=2*order;
  rl* glpts=new rl[glorder];
  rl* glwts=new rl[glorder];
  gltable(glorder,glpts,glwts);
  for(int glk=0;glk<glorder;glk++){
    rl* lvals=legendrevalues(order,glpts[glk]);
    for(int i=0;i<order;i++){
      for(int j=0;j<order;j++){
	retmat[arrayindx(i,j,order,order)]+=lvals[i]*lvals[j]*glwts[glk]*
	  func(ytox(glpts[glk],this,dbasis))*dxdy;
      }
    }
  }
  return retmat;
}

cmplx* legendrebasis::evalbasisfuncs(rl y){
  rl* legvals=legendrevalues(order-1,y);
  //cout << "legvals\n";
  //printmat(legvals,order,1);
  cmplx* retarray=initializearray(order);
  for(int i=0;i<order;i++){
    retarray[i]=legvals[i];
  }
  delete [] legvals;
  return retarray;		   
}

rl** legendrebasis::integrationtable(int glorder){
  //returns table of integration points, weights
  rl* glpts=new rl[glorder];
  rl* glwts=new rl[glorder];
  gltable(glorder,glpts,glwts);
  rl** rettable=new rl*[2];
  rettable[0]=glpts;
  rettable[1]=glwts;
  return rettable;
}

rl legendrebasis::mapy2range(rl y, rl x1, rl x2){
  rl xret=x1+((y+1.)/2.)*(x2-x1);
  return xret;
}

rl legendrebasis::mapfromchildrange(rl childx, rl childxmin, rl childxmax){
  rl yret=-1.+2.*(childx-childxmin)/(childxmax-childxmin);
  return yret;
}

rl legendrebasis::maptochildrange(rl y, rl childxmin, rl childxmax){
  rl childxret=childxmin+(childxmax-childxmin)*(y+1.)/2.;
  return childxret;
}

rl legendrebasis::dy2dx(rl x1, rl x2){
  return (x2-x1)/2.;
}

//cmplx* legendrebasis::xpsi(cmplx* __restrict psi){
//  //multiply psi by x.  In terms of the local coordinates, this is 
//  //a+b*y, where a=(xmin+xmax)/2 and b=(xmax-xmin)/2
//  //Here it is very useful to use Bonnet's recursion formula
//  //(2n+1)x P_n=(n+1)P_{n+1} + n P_{n-1}
//
//  cmplx* __restrict retarray=initializerray(order);
//  rl a=(xmin+xmax)/2.;
//  rl b=(xmax-xmin)/2.;
//  for(int n=0;n<order;n++){
//    retarray[n]+=psi[n]*b;
//  }
//  for(int n=1;n<order;n++){
//    retarray[n]+=psi[n-1]*a*rl(n)/(2.*rl(n)+1.);
//  }
//  for(int n=0;n<(order-1);n++){
//    retarray[n]+=psi[n+1]*a*(rl(n)+1.)/(2.*rl(n)+1.);
//  }
//  return retarray;
//}
