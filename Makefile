CC = clang++ # icpc #g++ #icpc

#CMPLFLAGS=-lgsl -lgslcblas -llapack -lblas -lg2c -mkl=sequential -p -g -pg -fno-inline # for debugging
#CMPLFLAGS = -O3  -lgsl -lgslcblas -llapack -lblas -lg2c -mkl=sequential -lpthread #-p -g -fno-inline #-L$(MKLROOT)/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread #-lm #-mkl 

#CMPLFLAGS=-O3 -lgsl -lgslcblas -llapack -lblas -lm -lg2c -lpthread -std=c++0x #for gcc

CMPLFLAGS=-std=c++11 -O3 -framework Accelerate 
#CMPLFLAGS=-O3 -lgsl -lgfortran -lgslcblas  -lblas -llapack -lfftw3 -mkl=parallel -lpthread -std=c++0x -vec-report1 #-pg # -p -g -pg
#CMPLFLAGS=-O3 -lgsl -lgfortran -lgslcblas  -lblas -llapack -lfftw3 -mkl=parallel  -lpthread -std=c++0x -vec-report1 -p -g -pg #for debugging


#LAPACK = -llapack #-L/usr/local/intel/mkl_10.1.0.009beta/lib/em64t -L/usr/lib64 -lmkl_core -liomp5 -lmkl_solver -lmkl_lapack -lmkl_em64t -lmkl_intel_thread -lpthread -Vaxlib  -L/usr/local/cernlib/lib/ -lmathlib -lkernlib
#MKLROOT=/usr/local/intel/mkl/10.0.011
CMPLOPTIONS=-I$(MKLROOT)/include
#MKLROOT=/usr/local/intel/mkl/10.0.011/mkl/lib/em64t/
MKLPATH= #/opt/intel/Compiler/11.1/064/mkl/lib/em64t/
MKLINCLUDE= #/opt/intel/Compiler/11.1/064/mkl/lib/
LAPACK =


FFTW= #/usr/local/fftw/lib/libfftw3.a
RKF45= #/home/zwalters/fortran_library/RKF45/*.o


#ARPACK= /usr/lib/libarpack.so.2.1 #/usr/local/lib/libarpack.a
#ARPACK=#/usr/local/lib/libarpack_LNX-64.a

SPLINES = # ~zwalters/splines/*.a	

DEL=*.o onedpropagator.x
#SOURCE=gauleg.f zPackfgmres.f 
SOURCE=onedpropagator.cpp potential.cpp pulse.cpp legendrebasis.cpp legendre.cpp laguerre.cpp gltable.cpp laguerretable.cpp printmat.cpp basischange.cpp derivedbasis.cpp temporalbasis.cpp triDbasis.cpp laguerrebasis.cpp mapintegrationtable.cpp ecsconversion.cpp globalbasis.cpp arrayfunctions.cpp basis.cpp linearsystemsetup.cpp linearsystemsetup_lagrange.cpp globallinearsystemsetup.cpp globallinearsystemsetup_lagrange.cpp globalbasis_minimizeaction.cpp globalbasis_psi0setup.cpp wf.cpp wf_functions.cpp functions.cpp wfxt_minimizeaction.cpp wfxt_propagate.cpp wf_iterativesolution.cpp wf_iterativesolution_incoherentsource.cpp wf_iterativesolution_diagonalU.cpp wf_iterativesolution_temporaleigenvectors.cpp wf_iterativesolution_temporaleigenvectors_aux.cpp globalbasis_eltbdysetup.cpp sourcewf_inhomogeneous_action.cpp  wf_boundstate.cpp globalbasis_boundstate.cpp #wf_fftw.cpp
HEADER=classes.h 
MODSOURCE=

OBJECT=$(SOURCE:.cpp=.o)
#OBJECTF90=$(SOURCEF90:.f90=.o)
#MODOBJECT=$(MODSOURCE:.f90=.o)
#MODMOD=$(MODSOURCE:.o=.mod)

.SUFFIXES: .cpp .h 
onedprop.x: $(OBJECT) $(HEADER)
	$(CC) $(CMPLFLAGS) -o onedpropagator.x  $(OBJECT) 

.cpp.o: $(HEADER)            
	$(CC) $(CMPLFLAGS) -c $< 
.h.o: 
	$(CC) $(CMPLFLAGS)  -c $<
clean:
	rm $(DEL) 
