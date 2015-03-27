#include "./classes.h"



rl laguerremap(rl x, rl R0){
  if(R0>0){
    return x-R0;
  }
  return R0-x;
}

rl* laguerrevalues(int order, rl x){
  if(order<0){
    return 0;
  }
  rl* __restrict vals=new rl[order+1];
  vals[0]=1.;
  vals[1]=1.-x;
  for(int i=2;i<order+1;i++){
    rl rk=i-1.;
    vals[i]=((2.*rk+1.-x)*vals[i-1]-rk*vals[i-2])/(rk+1.);
  }
  return vals;
}

rl* laguerrederivs(int order,int alpha, int derivorder, rl x){
  if(order<0){
    return 0;
  }
  rl* __restrict vals=new rl[order+1];
  for(int i=0;i<order+1;i++){
    vals[i]=laguerrederiv(i,alpha,derivorder,x);
  }
  return vals;
}


rl laguerreval(int order, rl x){
  //find value of laguerre polynomial using three term recursion relation
  if(order<0){
    return 0;
  }
  if(order==0){
    return 1.;
  }
  if(order==1){
    rl retval=1.-x;
    return retval;
  }
  rl vm1=laguerreval(order-1,x);
  rl vm2=laguerreval(order-2,x);
  rl rk=order-1.;
  return ((2.*rk+1.-x)*vm1-rk*vm2)/(rk+1.);
}

rl laguerreval(int order, int alpha, rl x){
  //find value of associated laguerre polynomial by using three term recursion
  //relation in alpha (alpha must be integer)
  //x L_{n}^{alpha+1}(x)=(n+alpha)L_{n-1}^{alpha}(x)-(n-x)L_{n}^{alpha}(x)
  if(alpha<0){
    return 0.;
  }
  if(alpha==0){
    return laguerreval(order,x);
  }
  
  rl val1=laguerreval(order-1,alpha-1,x);
  rl val2=laguerreval(order,alpha-1,x);
  rl retval=(1.*order+1.*(alpha-1.))*val1-(1.*order-x)*val2;
  retval=retval/x;
  return retval;

}

rl laguerrederiv(int order, int alpha, int derivorder,rl x){
  return pow(-1.,derivorder)*laguerreval(order-derivorder,alpha+derivorder,x);
}
