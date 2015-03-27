#include "./classes.h"
//set up common array functions as BLAS calls

cmplx* arraycopy(cmplx* inparray,int order){
  cmplx* retarray=new cmplx[order];
  cblas_zcopy(order,inparray,1,retarray,1);
  return retarray;
}

cmplx* arrayconj(cmplx* inparray,int order){
  cmplx* retarray=new cmplx[order];
  cblas_zcopy(order,inparray,1,retarray,1);
  for(int i=0;i<order;i++){
    retarray[i]=conj(retarray[i]);
  }
  return retarray;
}

void arrayscale(cmplx* inparray,int order,cmplx scale){
  cblas_zscal(order,&scale,inparray,1);
}

void zeroarray(cmplx* inparray,int order){
  //arrayscale(inparray,order,0.);
  for(int i=0;i<order;i++){
    inparray[i]=0.;
  }
}

cmplx* initializearray(int order){
  cmplx* retarray=new cmplx[order];
  zeroarray(retarray,order);
  return retarray;
}

cmplx* arraymultiply(cmplx* inparray, int order, cmplx val){
  cmplx* retarray=arraycopy(inparray,order);
  arrayscale(retarray,order,val);
  return retarray;
}

cmplx* transpose(cmplx* inparray,int nx1, int nx2){
  cmplx* retarray=new cmplx[nx1*nx2];
  for(int i=0;i<nx1;i++){
    for(int j=0;j<nx2;j++){
      retarray[arrayindx(j,i,nx2,nx1)]=inparray[arrayindx(i,j,nx1,nx2)];
    }
  }
  return retarray;
}

cmplx* arraysum(cmplx coeff1, cmplx* __restrict array1, cmplx coeff2, 
		cmplx* __restrict array2,int nx){
  cmplx* retarray=arraymultiply(array1,nx,coeff1);
  cblas_zaxpy(nx,&coeff2,array2,1,retarray,1);
  //for(int i=0;i<nx;i++){
  //  retarray[i]=coeff1*array1[i]+coeff2*array2[i];
  //}
  return retarray;
}

cmplx* arraydiff(cmplx coeff1, cmplx* __restrict array1, cmplx coeff2, 
		cmplx* __restrict array2,int nx){
  cmplx coeff2p=-coeff2;
  return arraysum(coeff1,array1,coeff2p,array2,nx);
}

cmplx* arraysum(cmplx* __restrict array1, cmplx* __restrict array2, int nx){
  cmplx coeff1=1.;
  cmplx coeff2=1.;
  return arraysum(coeff1,array1,coeff2,array2,nx);
}

cmplx* arraydiff(cmplx* __restrict array1, cmplx* __restrict array2, int nx){
  cmplx coeff1=1.;
  cmplx coeff2=-1.;
  return arraysum(coeff1,array1,coeff2,array2,nx);
}

cmplx* square_banded_matrix_multiply(int nx, int kl1, int ku1, 
				     cmplx* __restrict array1, 
				     int kl2,int ku2,cmplx* __restrict array2){
  //multiply two square banded matrices & return the result in banded format
  int kunew=ku1+ku2;
  int klnew=kl1+kl2;
  int widthnew=2*klnew+kunew+1;
  cmplx* __restrict retarray=initializearray(nx*widthnew);
  for(int i=0;i<nx;i++){
    int jmin=max(0,i-kl1);
    int jmax=min(i+ku1,nx-1);
    //cout << "jmin, jmax "<<jmin<<", "<<jmax<<"\n";
    for(int j=jmin;j<=jmax;j++){
//      cmplx Aij=array1[bandarrayindx(i,j,ku1,kl1,nx)];
//      int kstart=max(0,j-kl2);
//      int kstop=min(j+ku2,nx-1);
//      int tmplength=kstop-kstart+1;
//      int array2startindx=bandarrayindx(j,kstart,ku2,kl2,nx);
//      int startindx=bandarrayindx(i,kstart,kunew,klnew,nx);
//      int inc=1;
//      cblas_zaxpy(tmplength,&Aij,array2+array2startindx,inc,
//		  retarray+startindx,inc);
//
      for(int k=max(0,j-kl2);k<=min(j+ku2,nx-1);k++){
	//cout << "i,j,k "<<i<<", "<<j<<", "<<k<<"\n";
	retarray[bandarrayindx(i,k,kunew,klnew,nx)]+=
	  array1[bandarrayindx(i,j,ku1,kl1,nx)]*
	  array2[bandarrayindx(j,k,ku2,kl2,nx)];
      }
    }
  }
  return retarray;
  
}
