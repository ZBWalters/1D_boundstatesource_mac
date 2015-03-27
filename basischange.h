//#include "./basis.h"
//#include "./legendrebasis.h"
//#include "./derivedbasis.h"
//#include "./laguerrebasis.h"

//functions for mapping from one basis to another

//many of these functions are actually introduced in classes.h so as to be
//more easily available in different parts of the program
rl** mapintegrationtable(legendrebasis* bas1, derivedbasis* bas2,
			 rl** integrationtable, int order);
rl xtoy(rl x, legendrebasis* bas1,derivedbasis* bas2);
void ecsconversion(derivedbasis bs, rl theta);
cmplx* borderfunctionbasis(int order);
cmplx* leftborderfunctionbasis(int order);
cmplx* rightborderfunctionbasis(int order);
cmplx* leftborderfunctionbasis_laguerre(int order);
cmplx* rightborderfunctionbasis_laguerre(int order);
cmplx* leftbordervec(int ox);
cmplx* rightbordervec(int ox);
cmplx* Cnvec(int n, int ox);
cmplx* Cnvec2(int n, int ox);
cmplx* identitybasis(int order);
cmplx* matrix_basischange(int oldorder, int neworder, 
			  cmplx* matrix, cmplx* bfmat);
cmplx* vector_basischange(int oldorder, int neworder,cmplx* vector, 
			  cmplx* bfmat);
rl map2intervals(rl x, rl xmin, rl xmax, rl ymin, rl ymax);
rl displacex(rl x, rl R0);
