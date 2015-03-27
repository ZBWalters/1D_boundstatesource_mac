#include "./classes.h"
//#include "./globalbasis.h"
#include "./linearsystemsetup.h"
//#include "./temporalbasis.h"

//set up linear systems to be solved, using lagrange multipliers to specify
//initial conditions.

//cmplx*** globalbasis::linearsystemarrays_lagrange(cmplx* psi_global){
//  //given psi_global, returns arrays containing the least action linear
//  //systems for every element in the global basis.
//
//  //first set up linear systems in each subregion
//  cmplx** rhsvecs=new cmplx*[nelts];
//  cmplx** lhsmats=new cmplx*[nelts];
//
//  //set up psi0 in each basis
//
//  for(int i=0;i<nelts;i++){
//    cout << "i, nelts"<<i<<", "<<nelts<<"\n";
//    int tmporder=elementbases[i]->order;
//    cmplx* psi0=new cmplx[tmporder];
//    for(int j=0;j<tmporder;j++){
//      psi0[j]=psi_global[eltfirstindices[i]+j];
//    }
//    cout << "tbasis umat\n";
//    printmat(tbasis->umat,torder,torder);
//    cout << "tbasis qmat\n";
//    printmat(tbasis->qmat,torder,torder);
//    cmplx** linsys=linearsystemsetup_lagrange(elementbases[i],tbasis,psi0);
//    rhsvecs[i]=linsys[0];
//    lhsmats[i]=linsys[1];
//    delete [] psi0;
//    delete [] linsys;
//  }
//
//  cmplx*** retarray=new cmplx**[2];
//  retarray[0]=rhsvecs;
//  retarray[1]=lhsmats;
//  return retarray;
//
//}
cmplx*** globalbasis::linearsystemarrays_lagrange(cmplx* psi0){
  cmplx* psixt_global=psi0xF0(nbf,torder,torder,psi0);
  cmplx* dpsi0=initializearray(nbf);
  cmplx*** retarray=linearsystemarrays_lagrange(psixt_global,dpsi0);
  delete [] psixt_global;
  delete [] dpsi0;
  return retarray;
}
cmplx*** globalbasis::linearsystemarrays_lagrange(cmplx* psixt_global, 
						  cmplx* dpsi0){
  //given psi_global, returns arrays containing the least action linear
  //systems for every element in the global basis.
  
  //the convention used here is that psixt_global is a matrix with both
  //spatial and temporal basis coefficients, representing an initial guess at
  //the least action solution (this guess doesn't have to be good, since the
  //initial guess + calculated correction will be the same regardless of the
  //initial guess).  dpsi0 gives the initial condition for the _correction_
  //Thus, if we wanted to find the least action psixt corresponding to initial
  //guess psi0, we could do this in either of two ways:
  //1) set initial guess 0, dpsi0=psi0
  //2) initial guess = psi0(x)*F0(t), dpsi0=0


  cout << "linearsystemarrays_lagrange dpsi0\n";
  printmat(dpsi0,nbf,1);

  cout << "linearsystemarrays_lagrange psixt_global\n";
  printmat(psixt_global,nbf,torder);

  //first set up linear systems in each subregion
  cmplx** rhsvecs=new cmplx*[nelts];
  cmplx** lhsmats=new cmplx*[nelts];

  //set up psi0 in each basis

  for(int i=0;i<nelts;i++){


    int tmporder=elementbases[i]->order;
    cmplx* __restrict eltdpsi0=new cmplx[tmporder];
    //read dpsi0 into eltdpsi0
    for(int j=0;j<tmporder;j++){
      eltdpsi0[j]=dpsi0[eltfirstindices[i]+j];
    }
    
    //read psixt_global into eltpsixt
    cmplx* __restrict eltpsixt=new cmplx[tmporder*torder];
    //zeroarray(eltpsixt,tmporder*torder);
    for(int j=0;j<tmporder;j++){
      for(int m=0;m<torder;m++){
    	int xindx=eltfirstindices[i]+j;
    	eltpsixt[arrayindx(j,m,tmporder,torder)]=
    	  psixt_global[arrayindx(xindx,m,nbf,torder)];
	//cout << "nbf "<<nbf<<"\n";
	//cout << "j,m,psixt_global[xindx,m] "<<j<<", "<<m<<", "<<lsindx(xindx,m,nbf,torder)<<", "<<psixt_global[arrayindx(xindx,m,nbf,torder)]<<"\n";
      }
    }
  
    cout << "eltpsixt\n";
    printmat(eltpsixt,tmporder,torder);
    cout << "eltdpsi0\n";
    printmat(eltdpsi0,tmporder,1);


    //cout << "tbasis umat\n";
    //printmat(tbasis->umat,torder,torder);
    //cout << "tbasis qmat\n";
    //printmat(tbasis->qmat,torder,torder);
    cout << "linearsystemarrays_lagrange tmats #2\n";
    cout << "umat\n";
    printmat(tbasis->umat,torder,torder);
    cout << "qmat\n";
    printmat(tbasis->qmat,torder,torder);
    cout << "Amat\n";
    printmat(tbasis->Amat,torder,torder);
  
    

    cmplx** linsys=linearsystemsetup_lagrange(elementbases[i],tbasis,eltpsixt,eltdpsi0);
    cout << "linsys[0] "<<i<<"\n";
    printmat(linsys[0],torder,1);
    rhsvecs[i]=linsys[0];
    lhsmats[i]=linsys[1];
    delete [] eltdpsi0;
    delete [] eltpsixt;
    delete [] linsys;
  }

  cmplx*** retarray=new cmplx**[2];
  retarray[0]=rhsvecs;
  retarray[1]=lhsmats;
  return retarray;

}
cmplx** globalbasis::globallinearsystem_lagrange_banded(cmplx* psi_global){
  cmplx* psixt_global=psi0xF0(nbf,torder,torder,psi_global);
  cmplx* dpsi0=initializearray(nbf);
  cmplx** retarray=globallinearsystem_lagrange_banded(psixt_global,dpsi0);
  delete [] psixt_global;
  delete [] dpsi0;
  return retarray;
}


cmplx** globalbasis::globallinearsystem_lagrange_banded(cmplx* psixt_global,
							cmplx* dpsi0){
  //cout << "entering globallinearsystem_lagrange_banded\n";
  int torderp=torder+1;//number of equations per spatial element including
		       //lagrange multipliers
  int neqelt=xorder*(torderp);//number of equations per element
  int ku=neqelt;
  int kl=neqelt;
  int width=2*kl+ku+1;
  //cout << "neqelt,ku,kl,width "<<neqelt<<" "<< ku<<" "<<kl<<" "<<width<<"\n";
  int nxbases=nelts*xorder-(nelts-1);//number of x basis functions, less
				     //double counted boundary functions
  int neq=nxbases*(torderp);
  
  //allocate rhsvec, lhsmat for global linear system
  cmplx* __restrict globalrhsvec=initializearray(neq);
  cmplx* __restrict globallhsmat=initializearray(width*neq);

  //set up linear systems for each subelement
  cmplx*** eltlinsystems=linearsystemarrays_lagrange(psixt_global,dpsi0);
  cmplx** eltrhsvecs=eltlinsystems[0];
  cmplx** eltlhsmats=eltlinsystems[1];
  delete [] eltlinsystems;

  //copy element linear systems into global linear systems
  for(int eltnum=0;eltnum<nelts;eltnum++){
    //copy rhsvecs into globalrhsvec
    cout << "eltnum "<<eltnum<<"\n";
    cout << "eltrhsvec\n";
    printmat(eltrhsvecs[eltnum],xorder,torderp);
    for(int i=0;i<xorder;i++){
      int xindx=eltfirstindices[eltnum]+i;
      for(int n=0;n<torderp;n++){
	//cout << "global indx "<<i<<", "<<n<<", "<<lsindx(xindx,n,nxbases,torderp)<<"\n";
	//cout << "global value "<<globalrhsvec[lsindx(xindx,n,nxbases,torderp)]<<"\n";
	//cout << "local indx "<<i<<", "<<n<<", "<<lsindx(i,n,xorder,torderp)<<"\n";
	//cout << "local value "<<eltrhsvecs[eltnum][lsindx(i,n,xorder,torderp)]<<"\n";
	globalrhsvec[lsindx(xindx,n,nxbases,torderp)]+=
	  eltrhsvecs[eltnum][lsindx(i,n,xorder,torderp)];
      }
    }

    //cout << "rhsmat_banded\n";
    //printmat(globalrhsvec,neq,1);



    //    //copy lhsmats into globallhsmat
    //    for(int i=0;i<xorder;i++){
    //      for(int j=0;j<xorder;j++){
    //	for(int n=0;n<torderp;n++){
    //	  for(int m=0;m<torderp;m++){
    //	    int xindx1=eltfirstindices[eltnum]+i;
    //	    int xindx2=eltfirstindices[eltnum]+j;
    //	    int indx1=lsindx(xindx1,n,nxbases,torderp);
    //	    int indx1p=lsindx(i,n,xorder,torderp);
    //	    int indx2=lsindx(xindx2,m,nxbases,torderp);
    //	    int indx2p=lsindx(j,m,xorder,torderp);
    //	    globallhsmat[bandarrayindx(indx1,indx2,ku,kl,neq)]+=
    //	      eltlhsmats[eltnum][arrayindx(indx1p,indx2p,
    //					   xorder*(torderp),
    //					   xorder*(torderp))];
    //	  }
    //	}
    //      }
    //    }
    //  }

    //copy lhsmats into global lhsmat
    for(int i=0;i<xorder;i++){
      int xindx1=eltfirstindices[eltnum]+i;
      for(int j=0;j<xorder;j++){
	int xindx2=eltfirstindices[eltnum]+j;
	for(int n=0;n<torderp;n++){
	  for(int m=0;m<torderp;m++){
	    int indx1=lsindx(xindx1,n,nxbases,torderp);
	    int indx1p=lsindx(i,n,xorder,torderp);
	    int indx2=lsindx(xindx2,m,nxbases,torderp);
	    int indx2p=lsindx(j,m,xorder,torderp);
	    globallhsmat[bandarrayindx(indx1,indx2,ku,kl,neq)]+=
	      eltlhsmats[eltnum][arrayindx(indx1p,indx2p,
					   xorder*(torderp),
					   xorder*(torderp))];
	  }
	}
      }
    }
  }
	    

    cout<< "lhsmat_banded\n";
    printmat_banded(globallhsmat,neq,kl,ku,neq);

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
