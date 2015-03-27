#include "./classes.h"
#include "./basis.h"


//basis::basis(){
//  order=0;
//  omat=0;
//  nablamat=0;
//  nablasqmat=0;
//}
//
//basis::basis(int order_in){
//  order=order_in;
//  omat=new cmplx[order*order];
//  nablamat=new cmplx[order*order];
//  nablasqmat=new cmplx[order*order];
//  zeroarray(omat,order*order);
//  zeroarray(nablamat,order*order);
//  zeroarray(nablasqmat,order*order);
//}
//
//basis::basis(int order_in, cmplx* omat_in, cmplx* nablamat_in,
//	     cmplx* nablasqmat_in){
//  order=order_in;
//  omat=arraycopy(omat_in,order*order);
//  nablamat=arraycopy(nablamat_in,order*order);
//  nablasqmat=arraycopy(nablasqmat_in,order*order);
//}
//
//basis::basis(basis* inpbasis){
//  order=inpbasis->order;
//  omat=arraycopy(inpbasis->omat,order*order);
//  nablamat=arraycopy(inpbasis->nablamat,order*order);
//  nablasqmat=arraycopy(inpbasis->nablasqmat,order*order);
//}
//
//
//basis::~basis(){
//  delete [] omat;
//  delete [] nablamat;
//  delete [] nablasqmat;
//}

void basis::printomat(){
	if(omat!=0){
		printmat_scipy(omat,order,order);
	}
	else{
		cout << "null\n";
	}
}

void basis::printnablamat(){
	if(omat!=0){
		printmat_scipy(nablamat,order,order);
	}
	else{
		cout << "null\n";
	}
}

void basis::printnablasqmat(){
	if(omat!=0){
		printmat_scipy(nablasqmat,order,order);
	}
	else{
		cout << "null\n";
	}
}

cmplx* basis::oinvpsi(cmplx* psi){
	//find (omat)^-1 . psi
	cmplx* lhsmat=arraycopy(omat,order*order);
	cmplx* retvec=arraycopy(psi,order);

	int one=1;
	int* ipiv=new int[order];
	int info=0;
	zgesv_(&order, &one, reinterpret_cast <Lcmplx*> (lhsmat), &order, ipiv,
			reinterpret_cast <Lcmplx*> (retvec), &order, &info);
	if(info!=0){
		cout << "problem with zgesv_ in basis::oinvpsi!\n";
	}
	delete [] lhsmat;
	delete [] ipiv;

	return retvec;
}
