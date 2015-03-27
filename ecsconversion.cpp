#include "./classes.h"
//#include "./potential.h"
//#include "./basis.h"
//#include "./legendrebasis.h"
//#include "./laguerrebasis.h"
#include "./derivedbasis.h"

void ecsconversion(derivedbasis bs,rl theta){
  //given an exterior complex scaling angle theta, this subroutine adjusts the
  //matrices which make up the least action equation appropriately:
  //omat->omat*exp(I theta)
  //nablamat->nablamat
  //nablasqmat->nablasqmat*exp(-I theta)
  //Vmat->Vmat (if the long range potential goes as 1/r).

  cmplx II =cmplx (0.,1.);
  cmplx expitheta=exp(II*theta);
  cmplx expminitheta=exp(-II*theta);
  int ox=bs.order;
  for(int i=0;i<ox*ox;i++){
    bs.omat[i]*=expitheta;
    bs.nablasqmat[i]*=expminitheta;
  }
}

cmplx ecsX(rl x, rl recs, rl thetaecs){
  //if x>recs, return recs+(x-recs)*exp(II*thetaecs) 
  //if x<-recs, return -recs+(x+recs)*exp(II*thetaecs)
  cmplx II=cmplx(0.,1.);
  rl absR=abs(recs);
  cmplx retval=x;
  if(x>absR){
    retval=absR+(x-absR)*exp(II*thetaecs);
  }
  if(x<(-absR)){
    retval=-absR+(x+absR)*exp(II*thetaecs);
  }
  return retval;
}
