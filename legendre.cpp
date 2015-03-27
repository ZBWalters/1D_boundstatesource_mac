#include "./classes.h"

rl legmap(rl glx, rl t0, rl tf){
  //given glx in range [-1,1], map to t in [t0,tf]
  return t0+(glx+1.)*(tf-t0)/2.;
}

rl* legendrevalues(int order, rl x){
  // calculate values of legendre polynomials on range [-1,1]
  if(order<=0){
    return 0;
  }
  rl* __restrict vals=new rl[order+1];
  vals[0]=1.;
  if(order>=1){
    vals[1]=x;
  }
  for(int n=1;n<order;n++){
    //vals[n]=(2.*n+1.)*x*vals[n-1]-n*vals[n-2];
    vals[n+1]=((2.*n+1)*x*vals[n]-n*vals[n-1])/(n+1.);
  }
  return vals;
}
