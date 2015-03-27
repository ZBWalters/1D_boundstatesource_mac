//#include "./derivedbasis.h"
#include "./temporalbasis.h"

class globalbasis{
 public:
  int nelts;
  int nbdyelts;
  int xorder;//order of polynomials for x basis
  int torder;//order of polynomials for t basis
  rl laguerrekappa;
  rl ecstheta;
  rl ecsbdy;

  int nbf;
  
  cmplx* psi;
  cmplx* psi_xt;

  int* eltfirstindices;//gives the index of the global function basis
  //corresponding to the zeroth index of a particular function.
  
  legendrebasis* primitivetbasis;
  legendrebasis* primitivexbasis;
  laguerrebasis* primitivebdybasis_l;
  laguerrebasis* primitivebdybasis_r;
  

  temporalbasis* tbasis;//
  derivedbasis** elementbases;
  rl* elementbdys;
 

  

  globalbasis(int nelts_in, int nbdyelts_in,
	      rl dx,rl dt, int xorder_in,int torder_in, 
	      rl laguerrekappa_in, rl ecstheta_in);
  globalbasis(potential* pot,pulse* Pls, rl En_kappa, rl deltaphi, 
	      rl ecsbdy, rl ecsbuffer,rl dt, int xorder_in,int torder_in, 
	      rl laguerrekappa_in, rl ecstheta_in);
  ~globalbasis();
  void setupVmats(potential* pot);
  void updateHmats(potential* pot);
  void ecssetup(int nbdyelts_in, rl ecstheta_in);
  void ecssetup(int nbdyelts_in, rl ecstheta_in,potential* pot);

  int eltnum(rl x);

  rl nextdx(potential* pot,rl x, rl En, rl deltaphi);
  rl nextdx(potential* pot,rl x,rl ecsbdy, rl ecstheta,rl En, rl deltaphi);
  rl nextdx(potential* pot,pulse* Pls,rl x,rl ecsbdy, rl ecstheta,rl En, rl deltaphi);
  rl* elementboundaries_zerobased(potential* pot,pulse* Pls, 
				  rl En,rl deltaphi, rl ecsbdy, 
				  rl ecsbuffer,rl ecstheta,
				  int &nbdy,int &necselts, bool returnarray);
  rl* elementboundaries_zerobased(potential* pot,pulse* Pls,
				  rl En,rl deltaphi, rl ecsbdy, 
				  rl ecsbuffer, rl ecstheta, int &nbdy,int &necselts);

  //cmplx** globallinearsystemsetup();
  cmplx*** linearsystemarrays(cmplx* psi_global);
  cmplx*** linearsystemarrays_lagrange(cmplx* psixt_global, 
						    cmplx* dpsi0);
  cmplx** globallinearsystem_banded(cmplx* psiglobal);
  cmplx*** linearsystemarrays_lagrange(cmplx* psi_global);
  cmplx** globallinearsystem_lagrange_banded(cmplx* psiglobal);
  cmplx** globallinearsystem_lagrange_banded(cmplx* psixt_global,cmplx* dpsi0);

  cmplx* minimizeaction(cmplx* psi_global, rl t1, rl t2);
  cmplx* minimizeactioncorrection(cmplx* psixt_global, 
					       cmplx* dpsi0);
  cmplx* psi0setup(cmplx (*func)(rl x));
  cmplx* psi0setup(cmplx (*func)(rl x, rl kappa, rl rc),rl kappa, rl rc);
  cmplx* __restrict functimespsi(cmplx (*func)(rl x),cmplx* psi_in);

  rl norm(cmplx* psi_in);
  rl innerproduct(cmplx* psi1, cmplx* psi2);

  cmplx* eltLpsi(int eltnum);
  cmplx* eltLpsi(int eltnum,cmplx* psi_xt_in);

  cmplx* Olappsi(cmplx* psi_in);
  cmplx* Oinvpsi(cmplx* psi_in);
  cmplx* Hpsi(cmplx* psi_in);
  cmplx* ppsi(cmplx* psi_in);
  cmplx* Vpsi(cmplx* psi_in);
  cmplx* nablapsi(cmplx* psi_in);
  cmplx* nablasqpsi(cmplx* psi_in);
  cmplx* transposenablapsi(cmplx* psi_in);
  cmplx* globalLpsi();
  cmplx* globalLpsi(cmplx* psi_xt_in);
  cmplx psiLpsi(cmplx* psi_xt_1,cmplx* psi_xt_2);
  cmplx action();
  cmplx action(cmplx* psi_xt_in);

  cmplx* __restrict boundstate(int nstate, rl& En);
  
};


//useful generic function
cmplx* psi0xF0(int xorder, int torder, int torderp,  cmplx* psi0);
