//#include "mcheck.h"
#include <iostream>
#include <fstream>
using namespace std;
#include <complex>
#include <iostream>
#include <sstream>
#include "pthread.h"
#include <stdio.h>
//#include "math.h"
#include <math.h>
#define MKL_Complex16 std::complex<double>
//#include "gsl/gsl_cblas.h"
//#include "mkl.h"
//#include "mkl_lapack.h"

#include <vector>
#include <Accelerate/Accelerate.h>
//#include <vecLib/clapack.h>

//typedef int __CLPK_integer;
//typedef int __CLPK_logical;
//typedef float	 __CLPK_real;
//typedef double   __CLPK_doublereal;
//typedef struct { __CLPK_real r, i; } __CLPK_complex;
//typedef struct { __CLPK_doublereal r, i; } __CLPK_doublecomplex;
//typedef __CLPK_logical (*__CLPK_L_FP)();
//typedef  int __CLPK_ftnlen;

typedef double rl;
typedef complex <double> cmplx;

typedef __CLPK_doublecomplex Lcmplx;
typedef __CLPK_doublereal Lrl;

//typedef __CLPK_doublereal rl;
//typedef __CLPK_doublereal rl;
//typedef complex <double> cmplx;
//typedef complex<double> MKL_Complex16;
//typedef __CLPK_doublecomplex cmplx;
//typedef MKL_Complex16 cmplx;
typedef string str;
//typedef complex<double> cmplx;
//#define cmplx std::complex<double>;
//using namespace std;
//#include "pulse.h"
//#include "potential.h"
//#include "basis.h"
//#include "legendrebasis.h"
//#include "laguerrebasis.h"
//#include "derivedbasis.h"
//#include "./arrayindx.h"
#include "./arrayindx.h"

//useful constant
//rl Pi=3.14159265358979323846;

//functions for change of basis
cmplx* borderfunctionbasis(int order);
cmplx* leftborderfunctionbasis(int order);
cmplx* rightborderfunctionbasis(int order);
cmplx* leftborderfunctionbasis_laguerre(int order);
cmplx* rightborderfunctionbasis_laguerre(int order);
cmplx* leftbordervec(int ox);
cmplx* rightbordervec(int ox);
cmplx* Cnvec(int n, int ox);
cmplx* identitybasis(int order);
cmplx* matrix_basischange(int oldorder, int neworder, 
			  cmplx* matrix, cmplx* bfmat);
cmplx* vector_basischange(int oldorder, int neworder,cmplx* vector, 
			  cmplx* bfmat);

rl map2intervals(rl x, rl xmin, rl xmax, rl ymin, rl ymax);
rl displacex(rl x, rl R0);

////helper functions
int arrayindx(int i, int j, int ni, int nj);
int bandarrayindx(int i, int j, int ku, int kl, int nj);
int hbandarrayindx(char uplo, int i, int j, int ka);
int twodarrayindx(int i, int j, int n, int m, int ox, int ot);
int lsindx(int nx, int nt, int ox, int ot);
int lsindx2(int nx1, int nt1, int nx2, int nt2, int ox, int ot);



rl legmap(rl glx, rl t0, rl tf);
rl* legendrevalues(int order, rl x);
void gltable(int n, rl* glx, rl* glw);

rl laguerremap(rl x, rl R0);
rl* laguerrevalues(int order, rl x);
rl laguerreval(int order, rl x);
rl laguerreval(int order, int alpha, rl x);
rl laguerrederiv(int order, int alpha, int derivorder,rl x);
rl* laguerrederivs(int order,int alpha, int derivorder, rl x);
void laguerretable(int n, rl* glx, rl* glw);

void printmat(rl* mat, int n1, int n2);
void printmat(cmplx* mat, int n1, int n2);
void printmat_scipy(cmplx* mat, int n1, int n2);
void printmat_scipy(rl* mat, int n1, int n2);
void printmat_banded(cmplx* mat,int n1, int kl, int ku, int n2);
void printmat_formatted(str filename,cmplx* mat, int n1,int n2);
str matstr_formatted(cmplx* mat, int n1,int n2);
str matpairstr_formatted(rl* xmat, cmplx* zmat, int n);
void printmatpair(str filename,rl* xmat, cmplx* zmat, int n);

//array helper functions (wrappers for BLAS calls)
cmplx* arraycopy(cmplx* inparray,int order);
cmplx* arrayconj(cmplx* inparray,int order);
void arrayscale(cmplx* inparray,int order,cmplx scale);
void zeroarray(cmplx* inparray,int order);
cmplx* initializearray(int order);
cmplx* arraymultiply(cmplx* inparray, int order, cmplx val);
cmplx* transpose(cmplx* inparray,int nx1, int nx2);
cmplx* arraysum(cmplx coeff1, cmplx* __restrict array1, cmplx coeff2, 
		cmplx* __restrict array2,int nx);
cmplx* arraydiff(cmplx coeff1, cmplx* __restrict array1, cmplx coeff2, 
		 cmplx* __restrict array2,int nx);
cmplx* arraysum(cmplx* __restrict array1, cmplx* __restrict array2, int nx);
cmplx* arraydiff(cmplx* __restrict array1, cmplx* __restrict array2, int nx);
cmplx* square_banded_matrix_multiply(int nx, int kl1, int ku1, 
				     cmplx* __restrict array1, 
				     int kl2,int ku2,cmplx* __restrict array2);



//simple functions
rl factorial(int n);

cmplx gauss(rl x);
cmplx gauss(rl alpha, rl x);

cmplx ExpIwt(rl t);
cmplx ExpIwt(rl t, rl w0);

cmplx Expkappax(rl x, rl kappa, rl rc);
cmplx H1s(rl x, rl kappa, rl rc);
cmplx zH1s(rl x, rl kappa, rl rc);
cmplx X_func(rl x);

//exterior complex scaling functions
cmplx ecsX(rl x, rl recs, rl thetaecs);

//convert numeric type to string (template)
template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}


//matrix elements for tridiagonal basis
rl triD_U(int n, int m);
rl triD_U(int n, int m, rl dt);
rl triD_Q(int n, int m);
rl triD_Q(int n, int m, rl dt);

cmplx* triD_basischange(int order);
cmplx* triD_Umat(int order);
cmplx* triD_Qmat(int order);
cmplx* triD_Uinv(int order);
cmplx* triD_QUevals(int ot);
cmplx* triD_QUevecs(int ot);
cmplx** triD_QUeigensystem(int ot);
cmplx** generalized_eigensystem_reduced(int ot, cmplx* __restrict M1, 
					cmplx* __restrict M2);
cmplx* __restrict * __restrict padeigenvectors(int ot,cmplx* __restrict VR);

////functions for mapping from one basis to another
//rl** mapintegrationtable(legendrebasis* bas1, derivedbasis* bas2,
//			 rl** integrationtable, int order);
//rl xtoy(rl x, legendrebasis* bas1,derivedbasis* bas2);
//void ecsconversion(derivedbasis bs, rl theta);
