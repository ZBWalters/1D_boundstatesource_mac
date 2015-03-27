#include "./classes.h"
#include "./wf.h"


cmplx* __restrict * __restrict wf::inhomogeneous_action(rl En, temporalbasis* tbas){
  int ot=tbas->order;



  cmplx* __restrict Eiwtvec=tbas->Eiwtvec(En);//funcvector(ExpIwt,En);
  cout << "source En\t"<<En<<"\n";
  cout << "Eiwtvec\n";
  printmat(Eiwtvec,ot,1);

  cmplx* __restrict opsi=gbas->Oinvpsi(psi);
//  arrayscale(opsi,nbf,-1.);

  cmplx* __restrict * __restrict sourcepsixt=new cmplx*[ot];
  cmplx* __restrict * __restrict retaction=new cmplx*[ot];
  for(int n=0;n<ot;n++){
    retaction[n]=initializearray(nbf);
    cmplx* __restrict tmpsum=arraymultiply(psi,nbf,Eiwtvec[n]);
    //arrayscale(tmpsum,nbf,1.5);//pow(2.,.5));//test
    //arrayscale(tmpsum,nbf,-.5);
    arrayscale(tmpsum,nbf,-1.);
    sourcepsixt[n]=tmpsum;

  }
  for(int n=0;n<ot;n++){
    cmplx Aratio=0.;
    cmplx Asqratio=0.;
    deltapsiresidual_Aterms(n,ot,sourcepsixt[n],Aratio,Asqratio,retaction,tbas);
    //applydeltapsi(n,ot,sourcepsixt[n],sourcepsixt,retaction,tbas);
  }



  for(int n=0;n<ot;n++){
    delete [] sourcepsixt[n];
  }
  delete [] sourcepsixt;
  delete[] Eiwtvec;
  delete [] opsi;

  return retaction;
}

//cmplx* __restrict * __restrict wf::inhomogeneous_action(rl En, temporalbasis* tbas){
//  //return -Amat.pmat.wf+Asqmat.olap.wf -- ie, the inhomogenous source term
//  //for the least action equation
//  
//
//  int ot=tbas->order;
//  //calculate P.psi for A.p terms
//  cmplx* __restrict pwf=gbas->ppsi(psi);
//  cmplx* __restrict opwf=gbas->Olappsi(pwf);
//  //calculate Olap.psi for A^2 term
//  cmplx* __restrict olapwf=gbas->Olappsi(psi);
//  cmplx* __restrict olapsqwf=gbas->Olappsi(olapwf);
//  cmplx* __restrict powf=gbas->ppsi(olapwf);
//  //the time evolution of the bound state is given by Exp(-I E0 t)
//  //cmplx* __restrict Eiwtvec_func=tbas->funcvector(ExpIwt,-En);
//
//
//  
//
//  
//
//  //call zgemv to find Amat.ExpIwtvector
//  //Eiwtvec is integral of function against basis sets, so goes as deltat
//  //cmplx* __restrict Eiwtvec=initializearray(tbas->order);
//  cmplx* __restrict Eiwtvec=tbas->funcolap(ExpIwt,En);
//  cmplx* __restrict AdotEiwt=initializearray(tbas->order);
//  cmplx* __restrict AsqdotEiwt=initializearray(tbas->order);
//  cmplx alpha=1.;
//  cmplx beta=0.;
////  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
////	      &alpha,tbas->umat,ot,Eiwtvec_func,1,&beta,Eiwtvec,1);
//  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
//	      &alpha,tbas->Amat,ot,Eiwtvec,1,&beta,AdotEiwt,1);
//
////  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
////	      &alpha,tbas->Amat,ot,AdotEiwt,1,&beta,AsqdotEiwt,1);
//  
//  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
//	      &alpha,tbas->Asqmat,ot,Eiwtvec,1,&beta,AsqdotEiwt,1);
//
//  //cout <<"Eiwtvec\n";
//  //printmat(Eiwtvec,1,ot);
//  //cout <<"AdotEiwt\n";
//  //printmat(AdotEiwt,1,ot);
//  //cout <<"AsqdotEiwt";
//  //printmat(AsqdotEiwt,1,ot);
//
//  //return A.p.psi.e^(-I w0 t)-1/2 Asq.Olap.psi..e^(-I w0 t)
//  cmplx* __restrict * __restrict retactionarray=new cmplx* __restrict [ot];
//  for(int n=0;n<ot;n++){
//    //cmplx* __restrict tmparraypa=arraymultiply(powf,nbf,AdotEiwt[n]);
//    //cmplx* __restrict tmparrayap=arraymultiply(opwf,nbf,AdotEiwt[n]);
//    //cmplx* __restrict tmparrayaa=arraymultiply(olapsqwf,nbf,AsqdotEiwt[n]);
//    //cmplx* __restrict tmpsum=initializearray(nbf);
//    cmplx* __restrict tmparraypa=arraymultiply(pwf,nbf,AdotEiwt[n]);
//    cmplx* __restrict tmparrayap=arraymultiply(pwf,nbf,AdotEiwt[n]);
//    cmplx* __restrict tmparrayaa=arraymultiply(olapwf,nbf,AsqdotEiwt[n]);
//    cmplx* __restrict tmpsum=initializearray(nbf);
//    for(int i=0;i<nbf;i++){
//      tmpsum[i]=0.5*(tmparraypa[i]+tmparrayap[i]-tmparrayaa[i]);
//      //tmpsum[i]-=0.5*tmparray2[i];//add the two arrays
//    }
//    retactionarray[n]=tmpsum;
//    delete [] tmparraypa;
//    delete [] tmparrayap;
//    delete [] tmparrayaa;
//  }
//
//  delete [] pwf;
//  delete [] opwf;
//  delete [] powf;
//  delete [] olapwf;
//  delete [] olapsqwf;
//  //  delete [] Eiwtvec_func;
//  delete [] Eiwtvec;
//  delete [] AdotEiwt;
//  delete [] AsqdotEiwt;
//
//  return retactionarray;
//}


//////

cmplx* __restrict * __restrict wf
::inhomogeneous_action_lengthgauge(rl En, pulse* Pls, temporalbasis* tbas){
  return dipole_gauge(psi,En,Pls,tbas);
  //return timederiv_gauge(psi,En,Pls,tbas);
}

//cmplx* __restrict * __restrict wf::inhomogeneous_action_lengthgauge(rl En, pulse* Pls, temporalbasis* tbas){
//  //return E(t).zpsi, where input wf is assumed to be z*wf
//  
//
//  int ot=tbas->order;
//  rl dt=tbas->tmax-tbas->tmin;
//    //calculate Olap.psi 
//  cmplx* __restrict olapwf=gbas->Olappsi(psi);
//
//  cmplx* __restrict Eiwtvec=tbas->Eiwtvec(En);//funcvector(ExpIwt,En);//tbas->funcolap(ExpIwt,En);
//  //for(int n=0;n<ot;n++){
//  //  Eiwtvec[n]/=dt;
//  //}
//  cmplx* __restrict EdotEiwt=initializearray(tbas->order);//tbas->funcolap(ExpIwt,En);
//  
//  //call zgemv to find Emat.ExpIwtvector
//  //Eiwtvec is integral of function against basis sets, so goes as deltat
//  cmplx alpha=-1.;
//  cmplx beta=0.;
////  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
////	      &alpha,tbas->umat,ot,Eiwtvec_func,1,&beta,Eiwtvec,1);
////  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
////	      &alpha,tbas->Emat,ot,Eiwtvec,1,&beta,EdotEiwt,1);
//
//
//  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
//	      &alpha,tbas->Emat,ot,Eiwtvec,1,&beta,EdotEiwt,1);
//
//  cout << "Eiwtvec\n";
//  printmat(Eiwtvec,ot,1);
//  cout << "EdotEiwt\n";
//  printmat(EdotEiwt,ot,1);
//
//  //return E.zpsi.e^(-I w0 t)
//  cmplx* __restrict * __restrict retactionarray=new cmplx* __restrict [ot];
//  for(int n=0;n<ot;n++){
//    //cmplx* __restrict tmpsum=arraymultiply(olapwf,nbf,EdotEiwt[n]);
//    cmplx* __restrict tmpsum=arraymultiply(psi,nbf,Eiwtvec[n]);
//    //arrayscale(tmpsum,nbf,1.5);//pow(2.,.5));
//    retactionarray[n]=tmpsum;
//  }
//
//  cmplx* __restrict * __restrict retactionarray_gauge=
//    dipole_gauge(retactionarray,Pls,tbas);
//
//  delete [] olapwf;
//  delete [] Eiwtvec;
//  delete [] EdotEiwt;
//
//  for(int n=0;n<ot;n++){
//    delete [] retactionarray[n];
//  }
//  delete [] retactionarray;
//
//  return retactionarray_gauge;
//}

//cmplx* __restrict * __restrict wf::inhomogeneous_action_lengthgauge_gaugeterm(rl En, pulse* Pls,temporalbasis* tbas){
//  //calculates the inhomogeneous action due to a source term zpsi which is
//  //calculated in the length gauge.  Because of the change of gauge, there is
//  //an exponential term exp[i x A(x,t)] multiplying F cos(w t) zpsi exp(i E0 t)
//
//  cmplx II=cmplx(0.,1.);
//  
//
//  int torder=tbas->order;
//  rl** tinttable=tbas->integrationtable(torder);
//  rl* __restrict tintpts=tinttable[0];
//  rl* __restrict tintwts=tinttable[1];
//
//  //initialize return arrays
//  cmplx* __restrict * __restrict retarray=new cmplx*[torder];
//  for(int n=0;n<torder;n++){
//    retarray[n]=initializearray(nbf);
//  }
//
//  cmplx* __restrict * __restrict tbfvals=new cmplx*[torder];
//  for(int n=0;n<torder;n++){
//    tbfvals[n]=tbas->evalbasisfuncs(tintpts[n]);
//  }
//
//  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
//    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
//    int xorder=tmpbasis->order;
//    rl** xinttable=tmpbasis->integrationtable(xorder);
//    rl* __restrict xintpts=xinttable[0];
//    rl* __restrict xintwts=xinttable[1];
//
//    cmplx* __restrict * __restrict xbfvals=new cmplx*[xorder];
//    for(int i=0;i<xorder;i++){
//      xbfvals[i]=tmpbasis->evalbasisfuncs(tintpts[i]);
//    }
//
//    //integrate over temporal and spatial degrees of freedom
//    cmplx* __restrict integrals=initializearray(xorder*torder);
//    for(int i=0;i<xorder;i++){
//      for(int j=0;j<xorder;j++){
//	cmplx xpsi=tmpbasis->evalcoeffvector(psi+eltfirstindices[eltnum],
//					     xintpts[j]);
//	for(int n=0;n<torder;n++){
//	  for(int m=0;m<torder;m++){
//	    cmplx Aval=Pls->Az(tintpts[m]);
//	    cmplx Fval=Pls->Ez(tintpts[m]);
//	    integrals[arrayindx(i,n,xorder,torder)]+=
//	      conj(xbfvals[j][i]*tbfvals[m][n])*
//	      xpsi*Fval*exp(-II*En*tintpts[j])*//exp(-II*Aval*xintpts[j])*
//	      xintwts[j]*tintwts[m];
//	  }
//	}
//	
//
//      }
//    }
//    //copy integrals into retarray
//    for(int n=0;n<torder;n++){
//      for(int i=0;i<xorder;i++){
//	retarray[n][eltfirstindices[eltnum]+i]+=
//	  integrals[arrayindx(i,n,xorder,torder)];
//      }
//    }
//
//    delete [] xbfvals;
//    delete [] xintpts;
//    delete [] xintwts;
//    delete [] xinttable;
//
//    
//    
//  }
//    return retarray;
//}

////////////////////////////////////////////////
//cmplx* __restrict * __restrict wf::
//inhomogeneous_action_lengthgauge_dressed(rl En, pulse* Pls,temporalbasis* tbas){
//  //return E(t).zpsi, where input wf is assumed to be z*wf
//  
//
//  int ot=tbas->order;
//  rl dt=tbas->tmax-tbas->tmin;
//  rl tavg=(tbas->tmax+tbas->tmin)/2.;
//  rl Asqval=pow(Pls->Az(tavg),2.);
//    //calculate Olap.psi 
//  cmplx* __restrict olapwf=gbas->Olappsi(psi);
//  
//  rl dressedEn=En+Asqval;
//
//  cmplx* __restrict Eiwtvec=tbas->funcvector(ExpIwt,dressedEn);//tbas->funcolap(ExpIwt,En);
//  //for(int n=0;n<ot;n++){
//  //  Eiwtvec[n]/=dt;
//  //}
//  cmplx* __restrict EdotEiwt=initializearray(tbas->order);//tbas->funcolap(ExpIwt,En);
//  
//  //call zgemv to find Emat.ExpIwtvector
//  //Eiwtvec is integral of function against basis sets, so goes as deltat
//  cmplx alpha=-1.;
//  cmplx beta=0.;
////  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
////	      &alpha,tbas->umat,ot,Eiwtvec_func,1,&beta,Eiwtvec,1);
////  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
////	      &alpha,tbas->Emat,ot,Eiwtvec,1,&beta,EdotEiwt,1);
//
//
//  cblas_zgemv(CblasColMajor,CblasNoTrans,ot,ot,
//	      &alpha,tbas->Emat,ot,Eiwtvec,1,&beta,EdotEiwt,1);
//
//  cout << "Eiwtvec\n";
//  printmat(Eiwtvec,ot,1);
//  cout << "EdotEiwt\n";
//  printmat(EdotEiwt,ot,1);
//
//  //return E.zpsi.e^(-I w0 t)
//  cmplx* __restrict * __restrict retactionarray=new cmplx* __restrict [ot];
//  for(int n=0;n<ot;n++){
//    //cmplx* __restrict tmpsum=arraymultiply(olapwf,nbf,EdotEiwt[n]);
//    cmplx* __restrict tmpsum=arraymultiply(psi,nbf,EdotEiwt[n]);
//    
//    retactionarray[n]=tmpsum;
//  }
//
//  delete [] olapwf;
//  delete [] Eiwtvec;
//  delete [] EdotEiwt;
//
//  return retactionarray;
//}
