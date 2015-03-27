#include "./classes.h"
//#include "./globalbasis.h"
#include "./linearsystemsetup.h"
//#include "./temporalbasis.h"

cmplx*** globalbasis::linearsystemarrays(cmplx* psi_global){
  //given psi_global, returns arrays containing the least action linear
  //systems for every element in the global basis.

  //first set up linear systems in each subregion
  cmplx** rhsvecs=new cmplx*[nelts];
  cmplx** lhsmats=new cmplx*[nelts];

  //set up psi0 in each basis

  for(int i=0;i<nelts;i++){
    int tmporder=elementbases[i]->order;
    cmplx* psi0=new cmplx[tmporder];
    for(int j=0;j<tmporder;j++){
      psi0[j]=psi_global[eltfirstindices[i]+j];
    }
    cmplx** linsys=linearsystemsetup(elementbases[i],tbasis,psi0);
    rhsvecs[i]=linsys[0];
    lhsmats[i]=linsys[1];
    delete [] psi0;
    delete [] linsys;
  }

  cmplx*** retarray=new cmplx**[2];
  retarray[0]=rhsvecs;
  retarray[1]=lhsmats;
  return retarray;

}


cmplx** globalbasis::globallinearsystem_banded(cmplx* psi_global){
  int neqelt=xorder*(torder-1);//number of equations per element
  int ku=neqelt;
  int kl=neqelt;
  int width=2*kl+ku+1;
  int nxbases=nelts*xorder-(nelts-1);//number of x basis functions, less
				     //double counted boundary functions
  int neq=nxbases*(torder-1);
  
  //allocate rhsvec, lhsmat for global linear system
  cmplx* __restrict globalrhsvec=initializearray(neq);
  cmplx* __restrict globallhsmat=initializearray(width*neq);

  //set up linear systems for each subelement
  cmplx*** eltlinsystems=linearsystemarrays(psi_global);
  cmplx** eltrhsvecs=eltlinsystems[0];
  cmplx** eltlhsmats=eltlinsystems[1];
  delete [] eltlinsystems;

  //copy element linear systems into global linear systems
  for(int eltnum=0;eltnum<nelts;eltnum++){
    //copy rhsvecs into globalrhsvec
    for(int i=0;i<xorder;i++){
      int xindx=eltfirstindices[eltnum]+i;
      for(int n=0;n<torder-1;n++){
	globalrhsvec[lsindx(xindx,n,nxbases,torder-1)]+=
	  eltrhsvecs[eltnum][lsindx(i,n,xorder,torder-1)];
      }
    }


    //copy lhsmats into globallhsmat
    for(int i=0;i<xorder;i++){
      for(int j=0;j<xorder;j++){
	for(int n=0;n<torder-1;n++){
	  for(int m=0;m<torder-1;m++){
	    int xindx1=eltfirstindices[eltnum]+i;
	    int xindx2=eltfirstindices[eltnum]+j;
	    int indx1=lsindx(xindx1,n,nxbases,torder-1);
	    int indx1p=lsindx(i,n,xorder,torder-1);
	    int indx2=lsindx(xindx2,m,nxbases,torder-1);
	    int indx2p=lsindx(j,m,xorder,torder-1);
	    globallhsmat[bandarrayindx(indx1,indx2,ku,kl,neq)]+=
	      eltlhsmats[eltnum][arrayindx(indx1p,indx2p,
					   xorder*(torder-1),
					   xorder*(torder-1))];
	  }
	}
      }
    }
  }

  //delete rhsvecs, lhsmats corresponding to each element
  for(int i=0;i<nelts;i++){
    delete [] eltrhsvecs[i];
    delete [] eltlhsmats[i];
  }
  delete [] eltrhsvecs;
  delete [] eltlhsmats;

  cmplx** retarray=new cmplx*[2];
  retarray[0]=globalrhsvec;
  retarray[1]=globallhsmat;
  return retarray;
}
