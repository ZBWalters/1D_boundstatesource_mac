#include "./derivedbasis.h"
#include "./pulse.h"

class temporalbasis : public derivedbasis{
 public:
  cmplx*& umat;
  cmplx*& qmat;
  cmplx* Amat;
  cmplx* Asqmat;
  cmplx* Emat;

  rl& tmin;//This will be set to reference xmin of the parent basis
  rl& tmax;//This will be set to reference xmax of the parent basis

  //  temporalbasis(legendrebasis* prnt_in,rl dt);
  temporalbasis(legendrebasis* bas, cmplx* bfmat_in,rl tstart,rl tstop);
  temporalbasis(legendrebasis* bas, cmplx* bfmat_in,pulse* Pls,
		rl tstart,rl tstop);
  temporalbasis(temporalbasis* tbas, cmplx* bfmat_in,pulse* Pls,
		rl tstart,rl tstop);
  
  ~temporalbasis();

  

  void updateTmats(rl t1, rl t2, pulse* Pls);
  void Amatsetup(pulse* Pls);

  cmplx* Eiwtvec(rl En);
};