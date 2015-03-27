#include "./classes.h"
#include "./globalbasis.h"

cmplx* globalbasis::psi0setup(cmplx (*func)(rl x)){
	//given some initial function, find basis function coefficients A_{j} so
	// that \int_{-inf}^{inf} f(x) X_{i}(x)= 
	// \int_{-inf}^{inf} \sum_{j} A_{j} X_{j}(x)X_{i}(x)
	//or, in other words, \int_{-inf}^{inf} f(x) X_{i}(x)=\sum_{j}A_{j} O_{ij}
	//where O_{ij} is the overlap matrix


	//to accomplish this, first integrate func(x) against all of the primitive basis functions, then combine basis functions across element boundaries to set up the final linear system to be solved

	cmplx** __restrict primitiveoverlaps=new cmplx*[nelts];
	for(int eltnum=0;eltnum<nelts;eltnum++){
		int eltorder=elementbases[eltnum]->order;
		cmplx* __restrict eltoverlaps=new cmplx[eltorder];
		primitiveoverlaps[eltnum]=eltoverlaps;
		for(int i=0;i<eltorder;i++){
			eltoverlaps[i]=0.;
		}


		//integrationtable gives points and weights for gaussian quadrature
		rl** __restrict eltintegrationtable=elementbases[eltnum]->integrationtable(eltorder);
		rl* __restrict eltpts=eltintegrationtable[0];
		rl* __restrict eltwts=eltintegrationtable[1];

		for(int j=0;j<eltorder;j++){
			cmplx funcval=func(eltpts[j]);
			//cout << "funcval "<<eltpts[j]<<"\t"<<funcval<<"\n";
			cmplx* __restrict bfvals=elementbases[eltnum]->evalbasisfuncs(eltpts[j]);
			for(int i=0;i<eltorder;i++){
				eltoverlaps[i]+=bfvals[i]*funcval*eltwts[j];
			}
			delete [] bfvals;
		}

		//There's some numerical instability wrt boundary elements.  Try fixing this
		//by setting eltoverlaps to zero for these elements.
		//if((eltnum==0) or (eltnum==(nelts-1))){
		//  zeroarray(eltoverlaps,eltorder);
		//}

		//cout << "eltoverlaps \n";
		//printmat(eltoverlaps,eltorder,1);
		//for(int i=0;i<eltorder;i++){
		//  cout << eltoverlaps[i]<<"\t";
		//}
		//cout <<"\n";

		delete [] eltintegrationtable[0];
		delete [] eltintegrationtable[1];
		delete [] eltintegrationtable;


	}




	//Now combine functions in adjacent elements to create the total overlap
	//matrix in banded form used by zgbsv_

	int nxbases=nelts*xorder-(nelts-1);//number of x basis functions, less
	//double counted boundary functions
	int ku=xorder;
	int kl=xorder;
	int ldab=2*kl+ku+1;
	cmplx* __restrict bigomat=initializearray(nxbases*ldab);

	//now construct bigomat from the overlap matrices of the primitive basis
	//eltfirstindices corresponds to the index of the global function basis
	//which corresponds to a particular element's zeroth basis function
	for(int eltnum=0;eltnum<nelts;eltnum++){
		int eltorder=elementbases[eltnum]->order;
		int eltfirstindex=eltfirstindices[eltnum];
		cmplx* __restrict eltolap=elementbases[eltnum]->omat;
		//cout << "eltolap\n";
		//printmat(eltolap,eltorder,eltorder);
		for(int j=0;j<eltorder;j++){
			int globalj=eltfirstindex+j;
			for(int i=max(0,j-ku);i<min(eltorder,j+kl);i++){
				int globali=eltfirstindex+i;
				bigomat[bandarrayindx(globali,globalj,ku,kl,ldab)]+=
					//bigomat[arrayindx(kl+ku+globali-globalj,globalj,ldab,nxbases)]+=
					eltolap[arrayindx(i,j,eltorder,eltorder)];
			}
		}
	}
	//cout <<"bigomat\n";
	//printmat(bigomat,ldab,nxbases);


	//do the same thing to combine overlaps between the function and the primitive basis functions
	cmplx* __restrict retvec=initializearray(nbf);
	for(int eltnum=0;eltnum<nelts;eltnum++){
		int eltorder=elementbases[eltnum]->order;
		int eltfirstindex=eltfirstindices[eltnum];
		for(int i=0;i<eltorder;i++){
			retvec[eltfirstindex+i]+=primitiveoverlaps[eltnum][i];
		}
	}
	//cout << "retvec\n";
	//printmat(retvec,nbf,1);

	//now use the constructed overlap matrix and function overlaps to solve 
	// \int_{-inf}^{inf} f(x) X_{i}(x)=\sum_{j}A_{j} O_{ij}
	//where O_{ij} is the overlap matrix
	int neq=nbf;
	int nrhs=1;
	int* ipiv=new int[nbf];
	int info=0;
	zgbsv_(&neq, &kl, &ku, &nrhs, reinterpret_cast <Lcmplx*> (bigomat), &ldab,
			ipiv, reinterpret_cast <Lcmplx*> (retvec), &neq, &info);
	//cout << "zgbsv_ info "<<info<<"\n";
	if(info!=0){
		cout << "problem with zgbsv_ in psi0setup! "<<info<<"\n";
	}
	//cout << "solution retvec\n";
	//printmat(retvec,nbf,1);

	//delete primitiveoverlaps
	for(int eltnum=0;eltnum<nelts;eltnum++){
		delete [] primitiveoverlaps[eltnum];
	}
	delete [] primitiveoverlaps;
	delete [] bigomat;
	delete [] ipiv;


	return retvec;
}


cmplx* globalbasis::psi0setup(cmplx (*func)(rl x, rl kappa, rl rc),rl kappa, rl rc){
	//given some initial function, find basis function coefficients A_{j} so
	// that \int_{-inf}^{inf} f(x) X_{i}(x)= 
	// \int_{-inf}^{inf} \sum_{j} A_{j} X_{j}(x)X_{i}(x)
	//or, in other words, \int_{-inf}^{inf} f(x) X_{i}(x)=\sum_{j}A_{j} O_{ij}
	//where O_{ij} is the overlap matrix


	//to accomplish this, first integrate func(x) against all of the primitive basis functions, then combine basis functions across element boundaries to set up the final linear system to be solved

	cmplx** __restrict primitiveoverlaps=new cmplx*[nelts];
	for(int eltnum=0;eltnum<nelts;eltnum++){
		int eltorder=elementbases[eltnum]->order;
		cmplx* __restrict eltoverlaps=new cmplx[eltorder];
		primitiveoverlaps[eltnum]=eltoverlaps;
		for(int i=0;i<eltorder;i++){
			eltoverlaps[i]=0.;
		}


		//integrationtable gives points and weights for gaussian quadrature
		rl** __restrict eltintegrationtable=elementbases[eltnum]->integrationtable(eltorder);
		rl* __restrict eltpts=eltintegrationtable[0];
		rl* __restrict eltwts=eltintegrationtable[1];

		for(int j=0;j<eltorder;j++){
			cmplx funcval=func(eltpts[j],kappa,rc);
			//cout << "funcval "<<eltpts[j]<<"\t"<<funcval<<"\n";
			cmplx* __restrict bfvals=elementbases[eltnum]->evalbasisfuncs(eltpts[j]);
			for(int i=0;i<eltorder;i++){
				eltoverlaps[i]+=bfvals[i]*funcval*eltwts[j];
			}
			delete [] bfvals;
		}

		//There's some numerical instability wrt boundary elements.  Try fixing this
		//by setting eltoverlaps to zero for these elements.
		//if((eltnum==0) or (eltnum==(nelts-1))){
		//  zeroarray(eltoverlaps,eltorder);
		//}

		//cout << "eltoverlaps \n";
		//printmat(eltoverlaps,eltorder,1);
		//for(int i=0;i<eltorder;i++){
		//  cout << eltoverlaps[i]<<"\t";
		//}
		//cout <<"\n";

		delete [] eltintegrationtable[0];
		delete [] eltintegrationtable[1];
		delete [] eltintegrationtable;


	}




	//Now combine functions in adjacent elements to create the total overlap
	//matrix in banded form used by zgbsv_

	int nxbases=nelts*xorder-(nelts-1);//number of x basis functions, less
	//double counted boundary functions
	int ku=xorder;
	int kl=xorder;
	int ldab=2*kl+ku+1;
	cmplx* __restrict bigomat=initializearray(nxbases*ldab);

	//now construct bigomat from the overlap matrices of the primitive basis
	//eltfirstindices corresponds to the index of the global function basis
	//which corresponds to a particular element's zeroth basis function
	for(int eltnum=0;eltnum<nelts;eltnum++){
		int eltorder=elementbases[eltnum]->order;
		int eltfirstindex=eltfirstindices[eltnum];
		cmplx* __restrict eltolap=elementbases[eltnum]->omat;
		//cout << "eltolap\n";
		//printmat(eltolap,eltorder,eltorder);
		for(int j=0;j<eltorder;j++){
			int globalj=eltfirstindex+j;
			for(int i=max(0,j-ku);i<min(eltorder,j+kl);i++){
				int globali=eltfirstindex+i;
				bigomat[bandarrayindx(globali,globalj,ku,kl,ldab)]+=
					//bigomat[arrayindx(kl+ku+globali-globalj,globalj,ldab,nxbases)]+=
					eltolap[arrayindx(i,j,eltorder,eltorder)];
			}
		}
	}
	//cout <<"bigomat\n";
	//printmat(bigomat,ldab,nxbases);


	//do the same thing to combine overlaps between the function and the primitive basis functions
	cmplx* __restrict retvec=initializearray(nbf);
	for(int eltnum=0;eltnum<nelts;eltnum++){
		int eltorder=elementbases[eltnum]->order;
		int eltfirstindex=eltfirstindices[eltnum];
		for(int i=0;i<eltorder;i++){
			retvec[eltfirstindex+i]+=primitiveoverlaps[eltnum][i];
		}
	}
	//cout << "retvec\n";
	//printmat(retvec,nbf,1);

	//now use the constructed overlap matrix and function overlaps to solve 
	// \int_{-inf}^{inf} f(x) X_{i}(x)=\sum_{j}A_{j} O_{ij}
	//where O_{ij} is the overlap matrix
	int neq=nbf;
	int nrhs=1;
	int* ipiv=new int[nbf];
	int info=0;
	zgbsv_(&neq, &kl, &ku, &nrhs, reinterpret_cast <Lcmplx*> (bigomat), &ldab,
			ipiv, reinterpret_cast <Lcmplx*> (retvec), &neq, &info);
	//cout << "zgbsv_ info "<<info<<"\n";
	if(info!=0){
		cout << "problem with zgbsv_ in psi0setup! "<<info<<"\n";
	}
	//cout << "solution retvec\n";
	//printmat(retvec,nbf,1);

	//delete primitiveoverlaps
	for(int eltnum=0;eltnum<nelts;eltnum++){
		delete [] primitiveoverlaps[eltnum];
	}
	delete [] primitiveoverlaps;
	delete [] bigomat;
	delete [] ipiv;


	return retvec;
}



//  //given some initial function, find coefficients for the spatial basis
//  //functions so that psi0(x)=func(x)
//
//  //number of equations will be equal to number of basis functions (recalling
//  //that two adjacent elements share one basis function)
//  int neq=nelts*xorder+(nelts-1);//also equal to (nelts-1)*xorder+1
//
//  //set up arrays of function values 
//  cmplx** __restrict funcvals=new cmplx*[nelts];
//  cmplx** __restrict bfvals=new cmplx*[nelts];
//  for(int eltnum=0;eltnum<nelts;eltnum++){
//    int eltorder=elementbases[eltnum]->order;
//    rl** inttable=elementbases[eltnum]->integrationtable(eltorder);
//    rl* __restrict intpts=inttable[0];
//    rl* __restrict intwts=inttable[1];
//    cmplx* __restrict eltfuncvals=new cmplx[eltorder];
//    cmplx* __restrict eltbfvals=new cmplx[eltorder*eltorder];
//    //array of function values * intwts
//    for(int i=0;i<eltorder;i++){
//      eltfuncvals[i]=func(intpts[i]);
//    }
//    //array of basis function values * intwts
//    for(int i=0;i<eltorder;i++){
//      cmplx* tmpbfvals=elementbases[eltnum]->evalbasisfuncs(intpts[i]);
//      for(int j=0;j<eltorder;j++){
//	eltbfvals[arrayindx(i,j,eltorder,eltorder)]=tmpbfvals[j];
//      }
//      delete [] tmpbfvals;
//    }
//    funcvals[eltnum]=eltfuncvals;
//    bfvals[eltnum]=eltbfvals;
//
//    delete [] intpts;
//    delete [] intwts;
//    delete [] inttable;
//  }
//
//
//  //set up linear system to solve for coefficients & lagrange multipliers the
//  //convention is that the lagrange multiplier for the border functions comes
//  //in between the equation for the last element of the leftmost element and
//  //the zeroth equation for the rightmost element, so that the index for the
//  //first equation in a particular element is eltfirstindices[eltnum]+eltnum
//  int kl=xorder;
//  int ku=xorder;
//  int width=2*kl+ku+1;
//  cmplx* __restrict lhsmat=initializearray(width*neq);
//  cmplx* __restrict rhsvec=initializearray(neq);
//  for(int eltnum=0;eltnum<nelts;eltnum++){
//    int firstindx=eltfirstindices[eltnum]+eltnum;
//    int eltorder=elementbases[eltnum]->order;
//    //linear system:
//    //1)sum of coefficients * basis function values = function values
//    for(int i=0;i<eltorder;i++){
//      rhsvec[firstindx+i]=funcvals[eltnum][i];
//    }
//    for(int i=0;i<eltorder;i++){
//      for(int j=0;j<eltorder;j++){
//	lhsmat[bandarrayindx(firstindx+i,firstindx+j,ku,kl,neq)]=
//	  bfvals[eltnum][arrayindx(i,j,eltorder,eltorder)];
//      }
//    }
//  }
//  //set up equations for lagrange multipliers
//  for(int eltnum=0;eltnum<nelts-1;eltnum++){
//    int firstindx=eltfirstindices[eltnum]+eltnum;
//    int lagindx=firstindx+elementbases[eltnum]->order;
//    int leftindx=lagindx-1;//last indx in leftmost elt
//    int rightindx=lagindx+1;//zeroth indx in rightmost elt
//    lhsmat[bandarrayindx(leftindx,lagindx,ku,kl,neq)]+=1.;
//    lhsmat[bandarrayindx(lagindx,leftindx,ku,kl,neq)]+=1.;
//    lhsmat[bandarrayindx(rightindx,lagindx,ku,kl,neq)]-=1.;
//    lhsmat[bandarrayindx(lagindx,rightindx,ku,kl,neq)]-=1.;
//    rhsvec[lagindx]=0.;
//  }
//
//  //call zgbsv_ to solve for psi0
//  int nrhs=1;
//  int* ipiv=new int[neq];
//  int info=0;
//  zgbsv_(&neq,&kl,&ku,&nrhs,lhsmat,&width,ipiv,rhsvec,&neq,&info);
//  if(info!=0){
//    cout << "problem with zgbsv_ in psi0setup! "<<info<<"\n";
//  }
//
//  cmplx* retvec=initializearray(nbf);
//  for(int eltnum=0;eltnum<nelts;eltnum++){
//    int eltorder=elementbases[eltnum]->order;
//    for(int i=0;i<eltorder;i++){
//      retvec[eltfirstindices[eltnum]+i]=
//	rhsvec[eltfirstindices[eltnum]+eltnum+i];
//    }
//  }
//
//  for(int eltnum=0;eltnum<nelts;eltnum++){
//    delete [] funcvals[eltnum];
//    delete [] bfvals[eltnum];
//  }
//  delete [] funcvals;
//  delete [] bfvals;
//  
//  delete [] rhsvec;
//  delete [] lhsmat;
//  delete [] ipiv;
//  return retvec;
//
//
//}
