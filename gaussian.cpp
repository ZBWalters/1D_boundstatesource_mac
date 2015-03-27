#include "./classes.h"

cmplx gauss(rl x){
  cmplx II=cmplx(0.,1.);
  //cout << "gauss "<<x<<"\t"<< exp(II*x)*exp(-abs(pow(x,2)))<<"\n";
  //return exp(-abs(pow(x,2)));
  rl rc=.1;
  rl r=sqrt(pow(x,2)+pow(rc,2));
  //return exp(-abs(x));
  //return exp(-pow(x,2));
  return exp(1.*(-abs(r)+rc));
  
  //for rc large
  //return exp(-1.*(abs(r)-pow(rc,2)/(2*r)));

  //for rc small
  //return exp(-pow(x,2)/(2.*rc));
}

cmplx gauss(rl alpha, rl x){
  return exp(-alpha*pow(x,2));
}

cmplx ExpIwt(rl t){
  cmplx II=cmplx(0.,1.);
  rl w0=-.5;
  return exp(-II*w0*t);
}

cmplx ExpIwt(rl t, rl w0){
  cmplx II=cmplx(0.,1.);
  return exp(-II*w0*t);
}

cmplx Expkappax(rl x, rl kappa, rl rc){
  rl r=sqrt(pow(x,2)+pow(rc,2));
  return exp(kappa*r);
}

cmplx H1s(rl x, rl kappa, rl rc){
  rl r=sqrt(pow(x,2)+pow(rc,2));
  return exp(kappa*r);
}
