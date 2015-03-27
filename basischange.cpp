#include "./classes.h"
//#include "./basis.h"
//#include "./legendrebasis.h"
//#include "./laguerrebasis.h"
#include "./derivedbasis.h"
#include "./basischange.h"
//#include "./potential.h"

cmplx* leftbordervec(int ox){
  cmplx* __restrict evenvec=Cnvec2(0,ox);
  cmplx* __restrict oddvec=Cnvec2(1,ox);
  cmplx* __restrict psiret=new cmplx[ox];
  for(int i=0;i<ox;i++){
    psiret[i]=0.5*(evenvec[i]-oddvec[i]);//gives 1 at left bdy, 0 at right bdy
  }
  delete [] evenvec;
  delete [] oddvec;
  return psiret;
}

cmplx* rightbordervec(int ox){
  cmplx* __restrict evenvec=Cnvec2(0,ox);
  cmplx* __restrict oddvec=Cnvec2(1,ox);
  cmplx* __restrict psiret=new cmplx[ox];
  for(int i=0;i<ox;i++){
    psiret[i]=0.5*(evenvec[i]+oddvec[i]);//gives 0 at left bdy, 1 at right bdy
  }
  delete [] evenvec;
  delete [] oddvec;
  return psiret;
}


cmplx* Cnvec(int n, int ox){
  //returns a vector corresponding to the nth conforming basis function in the
  //x dimension, where a conforming basis function is one which is zero at the
  //left and right boundary.  The zeroth and first functions do not obey this
  //condition, but rather are constructed to be the even and odd functions
  //which are orthogonal to all conforming functions of a given order.

  cmplx* retvec=new cmplx[ox];
  for(int i=0;i<ox;i++){
    retvec[i]=0.;
  }
  if(n==0){//return even function orthogonal to all conforming functions
    cmplx sum=0.;
    for(int i=0;i<ox;i+=2){//even terms
      retvec[i]=(2.*i+1.);
      sum+=2.*i+1.;
    }
    //normalize such that Cnx[0](x)=1 at left and right boundaries
    for(int i=0;i<ox;i++){
      retvec[i]/=sum;
    }
    return retvec;
  }
  if(n==1){//return odd function orthogonal to all conforming functions
    cmplx sum=0.;
    for(int i=1;i<ox;i+=2){//odd terms
      retvec[i]=(2.*i+1.);
      sum+=2.*i+1.;
    }
    //normalize such that Cnx[1](x)=-1 at left boundary, 1 at right boundary
    for(int i=0;i<ox;i++){
      retvec[i]/=sum;
    }
    return retvec;
  }

      

  //try using even & odd border functions which aren't orthogonal to interior
  //functions.
  //if(n==0){
  //  retvec[0]=1.;
  //  return retvec;
  //}
  //if(n==1){
  //  retvec[1]=1.;
  //  return retvec;
  //}

  //for n!=0 or 1, Cn(x)=P_n(x)-P_0(x) for n even, P_n(x)-P_1(x) for n odd
  retvec[n]=1.;
  retvec[n%2]=-1.;
  return retvec;
}


cmplx* Cnvec2(int n, int ox){
  //second version of conforming basis functions.  Now even & odd boundary
  //functions are simply P_0 and P_1, while conforming basis functions are
  //given by BF_n(x)=.5(P_n(x)-P_{n-1}(x)), so that P_{2} and P_{3} have
  //nonzero overlap with the boundary functions.
  cmplx* retvec=initializearray(ox);
  if(n==0){
    retvec[0]=1.;
    return retvec;
  }
  if(n==1){
    retvec[1]=1.;
    return retvec;
  }
  //if neither of the above applies, return P_{n}-P_{n-2}
  retvec[n]=1.;
  retvec[n-2]=1.;
  return retvec;
}





cmplx* borderfunctionbasis(int order){
  //zeroth vector is left basis function
  //order-1st vector is right basis function
  //intermediate vectors give f_n(x)=P_n(x)-P_{mod(n,2)}(x)
  //for matrix basis[i,j], i is indx of basis function, 
  //j is indx of legendre function
  
  cmplx* retmat=initializearray(order*order);
  //cmplx* retmat=new cmplx[order*order];
  //zeroarray(retmat,order);
  
  //zeroth vector is left basis function
  int vecnum=0;
  //retmat[arrayindx(vecnum,0,order,order)]=.5;
  //retmat[arrayindx(vecnum,1,order,order)]=-.5;
  cmplx* __restrict lbdyvec=leftbordervec(order);
  for(int i=0;i<order;i++){
    retmat[arrayindx(vecnum,i,order,order)]=lbdyvec[i];
  }
  delete [] lbdyvec;

  //order-1st vector is left basis function
  vecnum=order-1;
  //retmat[arrayindx(vecnum,0,order,order)]=.5;
  //retmat[arrayindx(vecnum,1,order,order)]=.5;
  cmplx* __restrict rbdyvec=rightbordervec(order);
  for(int i=0;i<order;i++){
    retmat[arrayindx(vecnum,i,order,order)]=rbdyvec[i];
  }
  delete [] rbdyvec;
  
  for(int vecnum=1;vecnum<order-1;vecnum++){
    int rownum=vecnum+1;
    retmat[arrayindx(vecnum,rownum,order,order)]=1.;
    int colnum=rownum-2;
    //int colnum=rownum%2;
    retmat[arrayindx(vecnum,colnum,order,order)]=-1.;
  }
  
  //cout << "border function basis\n";
  //printmat(retmat,order,order);

  return retmat;
}

cmplx* leftborderfunctionbasis(int order){
  //for laguerre polynomials, let zeroth function have value 1 at left
  //boundary, all others have value 0.
  cmplx* retmat=initializearray(order*order);
  int vecnum=0;
  cmplx* __restrict lbdyvec=leftbordervec(order);
  for(int i=0;i<order;i++){
    retmat[arrayindx(vecnum,i,order,order)]=lbdyvec[i];
  }
  delete [] lbdyvec;
  //retmat[arrayindx(vecnum,0,order,order)]=1.;

  for(int vecnum=1;vecnum<order;vecnum++){
    retmat[arrayindx(vecnum,vecnum,order,order)]=1.;
    retmat[arrayindx(vecnum,0,order,order)]=-1.;
  }
  return retmat;
}



cmplx* rightborderfunctionbasis(int order){
  //for laguerre polynomials, let order-1st function have value 1 at right
  //boundary, all others have value 0.
  cmplx* retmat=initializearray(order*order);
  
  int vecnum=order-1;
  //retmat[arrayindx(vecnum,0,order,order)]=1.;
  cmplx* __restrict rbdyvec=rightbordervec(order);
  for(int i=0;i<order;i++){
    retmat[arrayindx(vecnum,i,order,order)]=rbdyvec[i];
  }
  delete [] rbdyvec;

  for(int vecnum=0;vecnum<order-1;vecnum++){
    retmat[arrayindx(vecnum,vecnum+1,order,order)]=1.;
    retmat[arrayindx(vecnum,0,order,order)]=-1.;
  }
  return retmat;
}

cmplx* leftborderfunctionbasis_laguerre(int order){
  //for laguerre polynomials, value at 0 is nCr(n+alpha,n), or 1 for alpha=0
  //choose linear combinations so that zeroth function has value 1 at 0, all
  //others have value 0 at 0.  I'll construct the functions with zero boundary
  //condition using sequential pairs of polynomials, so that the derived
  //overlap basis will be tridiagonal.
  cmplx* retmat=initializearray(order*order);
  //for 0 vector, let border function be 0th laguerre polynomial
  int vecnum=0;
  retmat[arrayindx(vecnum,0,order,order)]=1.;
  //for subsequent vectors, let border function be L_n(x)-L_{n-1}(x)
  for(int vecnum=1;vecnum<order;vecnum++){
    retmat[arrayindx(vecnum,vecnum,order,order)]=.5;
    retmat[arrayindx(vecnum,vecnum-1,order,order)]=-.5;
  }
  return retmat;
}
cmplx* rightborderfunctionbasis_laguerre(int order){
  //for laguerre polynomials, value at 0 is nCr(n+alpha,n), or 1 for alpha=0
  //choose linear combinations so that last function has value 1 at 0, all
  //others have value 0 at 0.  I'll construct the functions with zero boundary
  //condition using sequential pairs of polynomials, so that the derived
  //overlap basis will be tridiagonal.
  cmplx* retmat=initializearray(order*order);
  //for last vector, let border function be 0th laguerre polynomial
  int vecnum=order-1;
  retmat[arrayindx(vecnum,0,order,order)]=1.;
  //for prior vectors, let border function be L_n(x)-L_{n-1}(x)
  for(int i=1;i<order;i++){
    int vecnum=order-1-i;
    retmat[arrayindx(vecnum,i,order,order)]=.5;
    retmat[arrayindx(vecnum,i-1,order,order)]=-.5;
  }
  return retmat;
}



cmplx* identitybasis(int order){
  // returns identity matrix for the chosen order
  cmplx* retmat=initializearray(order*order);
  for(int i=0;i<order;i++){
    retmat[arrayindx(i,i,order,order)]=1.;
  }
  return retmat;
}


cmplx* matrix_basischange(int oldorder, int neworder, 
			  cmplx* matrix, cmplx* bfmat){
  //converts legmatrix calculated in original basis to matrix calculated wrt
  //new basis
  //Format of bfmat[i,j] is that i is index of new basis,
  //j is index of old basis

  cmplx* tmpmat=initializearray(oldorder*neworder);//new cmplx[oldorder*neworder];
  cmplx* retmat=initializearray(neworder*neworder);//new cmplx[neworder*neworder];

  cmplx alpha=1.;
  cmplx beta=0.;

  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasTrans,oldorder,neworder,oldorder,
	      &alpha,matrix,oldorder,bfmat,neworder,&beta,tmpmat,oldorder);
  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,neworder,neworder,oldorder,
	      &alpha,bfmat,neworder,tmpmat,oldorder,&beta,retmat,neworder);

//  char notrans='n';
//  char trans='t';
//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasTrans,order,order,order,
//	      &alpha,matrix,order,bfmat,order,&beta,tmpmat,order);
//  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,order,order,order,
//	      &alpha,bfmat,order,tmpmat,order,&beta,retmat,order);
  delete [] tmpmat;
  return retmat;
}



cmplx* vector_basischange(int oldorder, int neworder,cmplx* vector, 
			  cmplx* bfmat){
  cmplx* retvec=initializearray(neworder);
  cmplx alpha=1.;
  cmplx beta=0.;
  cblas_zgemv(CblasColMajor,CblasNoTrans,neworder,oldorder,&alpha,bfmat,
	      neworder,vector,1,&beta,retvec,1);
  return retvec;
}

rl map2intervals(rl x, rl xmin, rl xmax, rl ymin, rl ymax){
  return ymin+(x-xmin)/(xmax-xmin)*(ymax-ymin);
}

rl displacex(rl x, rl R0){
  return R0+x;
}

//rl ytox(rl y, legendrebasis lbas){
//  // map parent basis to derived basis
//  rl ymin=-1.;
//  rl ymax=1.;
//  
//  rl xret=xmin+((y-ymin)/(ymax-ymin))*(xmax-xmin);
//  return xret;
//}
//
//rl ytox(rl y, derivedbasis dbas1){
//  rl ymin=dbas1->xmin;
//  rl ymax=dbas1->xmax;
//
//  rl xret=xmin+((y-ymin)/(ymax-ymin))*(xmax-xmin);
//  return xret;
//}
//
//rl ytox(rl y, laguerrebasis lagbas){
//  
//}

