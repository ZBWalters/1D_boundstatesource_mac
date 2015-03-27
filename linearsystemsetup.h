//#include "./derivedbasis.h"
//#include "./temporalbasis.h"
#include "./globalbasis.h"



cmplx** tmatsetup_ics(int ot,cmplx* tmat, cmplx* umat,cmplx* icbfmat);
cmplx* Xpsi(int ix, int ox, int ot, cmplx* Xmat, cmplx* psi, cmplx coeff);
cmplx* Tpsi(int it, int ot, int ox, cmplx* Tmat, cmplx* psi, cmplx coeff);
cmplx* XTpsi(int ix, int ox, int it, int ot, cmplx* Xmat, cmplx* Tmat, 
	     cmplx* psi, cmplx coeff);
cmplx* Lpsi(int ox, int ot, cmplx** xmats, cmplx** tmats, 
	    cmplx* psi);
cmplx* XTmat(int ox1, int ox2, int ot1, int ot2, 
	     cmplx* Xmat, cmplx* Tmat, cmplx coeff);



cmplx* rhsvecsetup(int ot, int ox,cmplx** tmats, cmplx** xmats, cmplx* psi0);
cmplx* lhsmatsetup(int ot, int ox,cmplx** tmats, cmplx** xmats);
cmplx** linearsystemsetup(derivedbasis* xbas,temporalbasis* tbas,cmplx* psi0);

cmplx** linearsystemsetup_lagrange(derivedbasis* xbas,temporalbasis* tbas,
				   cmplx* psi0);
cmplx** linearsystemsetup_lagrange(derivedbasis* xbas,temporalbasis* tbas,
				   cmplx* psixt_guess,cmplx* dpsi0);
cmplx* lhsmatsetup_lagrange(int ot, int ox, cmplx** tmats, cmplx** xmats);
cmplx* rhsvecsetup_lagrange(int ot, int ox,cmplx* resid, cmplx* dpsi0);
cmplx* rhsvecsetup_lagrange(int ot, int ox, cmplx* resid);
//cmplx* psi0xF0(int xorder, int torder, int torderp,  cmplx* psi0);
