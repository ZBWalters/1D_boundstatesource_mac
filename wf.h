#include "./globalbasis.h"
//a wavefunction consists of an array of coefficients and a global basis to
//which they refer.
//The convention is that a wf contains only spatial coefficients, while
//a wfxt contains both spatial and temporal coefficients
class wf{
 public:
  int nbf;
  globalbasis* gbas;
  cmplx* psi;
  int* eltfirstindices;
  wf();
  wf(globalbasis* gbas_in);
  wf(globalbasis* gbas_in, cmplx* psi_in);
  wf(globalbasis* gbas_in, cmplx (*func)(rl x));
  ~wf();

  void propagate(pulse* Pls,rl tstart, rl tstop, rl dtstart,rl accuracygoal);

  void psi0setup(cmplx (*func)(rl x));
  void psi0setup(cmplx (*func)(rl x, rl kappa, rl rc),rl kappa, rl rc);
  void print();
  void print(int neltpts);
  str printstr(int neltpts);
  str printstr();

  rl psinorm();

  //methods for iterative method for propagation
  cmplx* __restrict * __restrict * __restrict  triD_iterativesolution(cmplx* __restrict  psi0, int bfordergoal, int dtordergoal,
			 int nmax, rl dt);
//  void correctbfsequence(int nmin, int nmax, cmplx** guessxt, 
//			 cmplx** residxt, rl dt);
  void correctbfsequence(int nmin, int nmax, cmplx* __restrict * __restrict  guessxt, 
			 cmplx* __restrict * __restrict  residxt, rl dt);
//  void correctioncycle(int n, cmplx** guessxt, cmplx** residxt, rl dt);
  void correctioncycle(int n, cmplx* __restrict * __restrict  guessxt,
		       cmplx* __restrict * __restrict  residxt, rl dt);
  void correctioncycle(int n,int nmax, cmplx* __restrict * __restrict  guessxt,
			 cmplx* __restrict * __restrict  residxt,
			 cmplx* __restrict * __restrict  targetxt,
		       temporalbasis* tbas);
  void correction_leveln(int n,int nmax,cmplx* __restrict * __restrict  guessxt,
			 cmplx* __restrict * __restrict  residxt,
			 cmplx* __restrict * __restrict  targetxt,
			 temporalbasis* tbas);
  void correct_differentlevels(int ncorr, int nresid, int nmax,
				 cmplx* __restrict * __restrict  guessxt,
				 cmplx* __restrict * __restrict  residxt,
				 cmplx* __restrict * __restrict  targetxt,
			       temporalbasis* tbas);
  void correctioncycle_dual_recursive(int nmax,int recursionlevel, 
					int maxrecursionlevel,rl accuracygoal,
					cmplx* __restrict * __restrict guessxt,
					cmplx* __restrict * __restrict residxt,
					cmplx* __restrict * __restrict targetxt,
				      temporalbasis* tbas);
  void correctionsequence_dual(int n1, int n2,int nmax, 
			       cmplx* __restrict * __restrict  guessxt,
			       cmplx* __restrict * __restrict  residxt,
			       cmplx* __restrict * __restrict  targetxt,
			       temporalbasis* tbas);
  void correctionsequence(int n1, int n2,int nmax, 
			  cmplx* __restrict * __restrict  guessxt,
			  cmplx* __restrict * __restrict  residxt,
			  cmplx* __restrict * __restrict  targetxt,
			  temporalbasis* tbas);

  cmplx* __restrict boundstate(int n, rl& En);
  void boundstatesetup(int nstate, rl& En);
  void boundstatesetup(int nstate, rl& En, rl t);
  void boundstatesetup_source(int nstate, rl& En);


  rl errornorm(int nmax,cmplx* __restrict * __restrict  residxt,
	       cmplx* __restrict * __restrict  targetxt);
    //  cmplx* residualcorrection(int n, cmplx** guesst, cmplx** residxt, rl dt);
  cmplx* residualcorrection(int n, cmplx* __restrict * __restrict  guesst, cmplx* __restrict * __restrict  residxt, rl dt);
  cmplx* residualcorrection(int n, cmplx* __restrict resid, 
			    temporalbasis* tbas);
  cmplx* residualcorrection(int ncorr,int nresid, cmplx* __restrict resid, 
		     temporalbasis* tbas);
//  void triD_applyupdate(int n, cmplx* deltapsin, 
//			cmplx** guessxt, cmplx** residxt, rl dt);
  void triD_updateresidual(int nmax,cmplx* __restrict * __restrict guessxt, 
			       cmplx* __restrict * __restrict residxt, rl dt);
  void triD_applydeltapsi(int n, cmplx* __restrict deltapsin, 
			cmplx* __restrict * __restrict guessxt, 
			cmplx* __restrict * __restrict residxt, rl dt);
  void applydeltapsi(int nmax,cmplx* __restrict * __restrict deltapsi,
			 cmplx* __restrict * __restrict guessxt,
			 cmplx* __restrict * __restrict residxt, 
		     temporalbasis* tbas );
  void applydeltapsi(int ncorr, int nmax, cmplx* __restrict deltapsi,
			 cmplx* __restrict * __restrict guessxt,
			 cmplx* __restrict * __restrict residxt, 
		     temporalbasis* tbas );
  void deltapsiresidual_Aterms(int ncorr,int nmax, cmplx* __restrict deltapsi,
				 cmplx Aratio, cmplx Asqratio,
			       cmplx* __restrict * __restrict residxt, 
			       temporalbasis* tbas );
  //  void printiterativesolution(int nmax, cmplx** guessxt, cmplx** residxt);
  void printiterativesolution(int nmax, cmplx* __restrict * __restrict  guessxt, cmplx* __restrict * __restrict  residxt);
  //  cmplx iterativeaction(cmplx** guessxt,cmplx** residxt, int nmax);
  cmplx iterativeaction(cmplx* __restrict * __restrict  guessxt,cmplx* __restrict * __restrict  residxt, int nmax);
  cmplx iterativeaction_inhomogenous(int nmax,
				       cmplx* __restrict * __restrict  guessxt,
				       cmplx* __restrict * __restrict  residxt,
				     cmplx* __restrict * __restrict  targetxt);
  void timestep(rl tmin,rl tmax,int steporder,rl accuracygoal,rl& newdt);
  //  cmplx* psif(cmplx** guessxt, int nmax);
  cmplx* psif(cmplx* __restrict * __restrict  guessxt, int nmax);
  //  cmplx* lastcomponentptr(cmplx** xtarray, int nmax);
  cmplx* lastcomponentptr(cmplx* __restrict * __restrict  xtarray, int nmax);
  //  rl lastcomponentsignificance(cmplx** xtarray, int nmax);
  rl lastcomponentsignificance(cmplx* __restrict * __restrict  xtarray, int nmax);
  rl nthcomponentsignificance(cmplx* __restrict * __restrict xtarray, int n);
  rl nthcomponentsignificance(cmplx* __restrict * __restrict xtarray1,
			      cmplx* __restrict * __restrict xtarray2, 
			      int n);
  cmplx* __restrict * __restrict 
    triD_legendre_conversion(int ot,cmplx* __restrict * __restrict tridwf);
  rl normsum_legendre(int ot,rl dt, 
		      cmplx* __restrict * __restrict legendrewf);
  rl lastlegendrecomponentsignificance(int ot,rl dt,
				       cmplx* __restrict * __restrict tridwf);
  void propagate_iterative(pulse* Pls, potential* Pot,int steporder,rl tmin,rl tmax, rl dtfirst, rl accuracygoal);



  //functions useful for propagation with inhomogenous source term
  cmplx* __restrict * __restrict inhomogeneous_action(rl En,
						      temporalbasis* tbas);
  cmplx* __restrict * __restrict 
    inhomogeneous_action_lengthgauge(rl En,pulse* Pls,temporalbasis* tbas);
  cmplx* __restrict * __restrict 
    inhomogeneous_action_lengthgauge_dressed(rl En,pulse* Pls, 
					     temporalbasis* tbas);
  cmplx* __restrict * __restrict 
    inhomogeneous_action_lengthgauge_gaugeterm(rl En, pulse* Pls,
					       temporalbasis* tbas);
//  cmplx* __restrict *__restrict 
//    dipole_gauge(cmplx* __restrict *__restrict wfxt, pulse* Pls,
//			temporalbasis* tbas);
  cmplx* __restrict *__restrict 
    dipole_gauge(cmplx* __restrict psi_in, rl En, pulse* Pls,
		 temporalbasis* tbas);
  cmplx* __restrict *__restrict 
    timederiv_gauge(cmplx* __restrict psi_in, rl En, pulse* Pls,
		    temporalbasis* tbas);
  cmplx* __restrict wf_gaugechange(cmplx* __restrict psi_in,
					pulse* Pls,rl t_in,rl plusminus);
  void psi_gaugechange(pulse* Pls,rl t_in,rl plusminus);
  void propagate_iterative_inhomogenous(rl tmin, rl tmax, rl& dt,
					pulse* Pls, potential* Pot,
					wf* sourcewf, rl sourceEn, 
					rl accuracygoal);
  void timestep(rl tmin, rl tmax, potential* Pot,pulse* Pls,
		wf* sourcewf, rl sourceEn, rl accuracygoal, rl& newdt);
  cmplx* __restrict * __restrict 
    residualdiff(cmplx* __restrict * __restrict guessxt,
		 cmplx* __restrict * __restrict targetxt,temporalbasis* tbas);

  //functions for use in Uinv iterative solution
  void initialguess_Uinv(int nmax, cmplx* __restrict psi0, 
		    cmplx* __restrict * __restrict psixt,
		    cmplx* __restrict * __restrict resid,
		    temporalbasis* tbas);
  void initialguess_Uinv(int nmax, cmplx* __restrict psi0, 
		    cmplx* __restrict * __restrict psixt,
		    cmplx* __restrict * __restrict resid,rl scale,
		    temporalbasis* tbas);
  cmplx* __restrict * __restrict Uinvpsi(int nmax,
					 cmplx* __restrict * __restrict dphi);
  cmplx* residualcorrection_Uinv(int ncorr,int nresid, cmplx* __restrict resid, 
			  temporalbasis* tbas);
  void applydeltapsi_Uinv(int nmax,cmplx* __restrict * __restrict deltapsi,
			 cmplx* __restrict * __restrict guessxt,
			 cmplx* __restrict * __restrict residxt, 
			  temporalbasis* tbas);
  void applydeltapsi_Uinv(int ncorr,int nmax, cmplx* __restrict deltapsi,
			 cmplx* __restrict * __restrict guessxt,
			 cmplx* __restrict * __restrict residxt, 
			  temporalbasis* tbas);
  void applydeltapsi_Uinv(int nmax,cmplx* __restrict * __restrict deltapsi,
			 cmplx* __restrict * __restrict guessxt,
			  cmplx* __restrict * __restrict residxt, rl scale,
			  temporalbasis* tbas);
  void applydeltapsi_Uinv(int ncorr,int nmax, cmplx* __restrict deltapsi,
			 cmplx* __restrict * __restrict guessxt,
			  cmplx* __restrict * __restrict residxt, rl scale,
			  temporalbasis* tbas);
  cmplx* __restrict * __restrict 
    cancelresiduals_Uinv(cmplx* __restrict * __restrict residxt,
			 temporalbasis* tbas);
  void timestep_Uinv(rl tmin, rl tmax, potential* Pot,pulse* Pls,
		wf* sourcewf, rl sourceEn, rl accuracygoal, rl& newdt);
  cmplx* __restrict * __restrict 
    psixtdiff(cmplx* __restrict * __restrict psi1, 
	      cmplx* __restrict * __restrict psi2, int ot);
  cmplx* __restrict * __restrict 
    deltaxt_Uinv(int recursionlvl, int maxrecursionlvl, rl accuracygoal,
		 cmplx * __restrict psi0,
		 cmplx* __restrict * __restrict targetxt,temporalbasis* tbas);
  rl maxdelta(cmplx* __restrict * __restrict deltaxt, int nmax);

  rl scalefactor(int ot, rl dt);
  void scalewfxt(cmplx* __restrict *__restrict wfxt,int ot, rl dt);
  void scalewfxt_inverse(cmplx* __restrict *__restrict wfxt,int ot, rl dt);
  cmplx* __restrict * __restrict  
    Qinvpsi(int nmax,cmplx* __restrict * __restrict dphi, rl scale);
  cmplx* __restrict * __restrict  
    Upsi(int nmax,cmplx* __restrict * __restrict dphi, rl scale);
  cmplx* __restrict * __restrict  
    Upsi(int nmax,cmplx* __restrict * __restrict dphi);
  cmplx* __restrict * __restrict psi_swaporders(int nmax,int order1,int order2,
		   cmplx* __restrict * __restrict dphi);
  cmplx* Olapinvpsi(int ncorr, cmplx* __restrict resid);
  void Upsi0(int nmax, cmplx* __restrict psi0, 
		    cmplx* __restrict * __restrict psixt,
		    cmplx* __restrict * __restrict resid,rl scale,
		    temporalbasis* tbas);

  cmplx* __restrict * __restrict 
    IOUinvQsolve(int ot,cmplx* __restrict * __restrict resid, rl scale);

  //functions for use in temporal eigenvector solve
  cmplx* __restrict deltawf_eigenvector(cmplx* __restrict resid, 
					cmplx Qval,cmplx Uval,
					cmplx Aval,cmplx Asqval);
  cmplx* __restrict deltawf_eigenvector_Olap_preconditioner(
					cmplx* __restrict resid, 
					cmplx Qval,cmplx Uval,
					cmplx Aval,cmplx Asqval);
  cmplx* __restrict commoneigenvectorcorrection_level0residuals(int nbf,int ot,
		         cmplx* __restrict resid0,
			 cmplx* __restrict * __restrict VR,cmplx Uval, 
			 cmplx Aval, cmplx Asqval);
  cmplx* __restrict * __restrict cancel0residuals_temporaleigenvector(
	 int nbf,int ot,
	 cmplx* __restrict resid0,
	 cmplx* __restrict * __restrict VR,
	 cmplx Uval,cmplx Aval, cmplx Asqval);
  cmplx* __restrict * __restrict canceleigenvectorresiduals(
	int nbf,int ot,cmplx* __restrict evals,
	cmplx* __restrict * __restrict padrevecs,
	cmplx* __restrict * __restrict padlevecs,
	cmplx* __restrict * __restrict targetxt,temporalbasis* tbas);
  cmplx* __restrict * __restrict canceleigenvectorresiduals_Olap_preconditioner(
       int nbf,int ot,cmplx* __restrict evals,
       cmplx* __restrict * __restrict padrevecs,
       cmplx* __restrict * __restrict padlevecs,
       cmplx* __restrict * __restrict targetxt,temporalbasis* tbas);
  cmplx* __restrict * __restrict eigenvectortowf(int nbf,int ot,
		cmplx* __restrict * __restrict dpsi,
		cmplx* __restrict * __restrict padrevecs);
  void applydeltapsi_temporaleigenvectors(int nbf,int nmax,
				   cmplx* __restrict * __restrict dpsi_ev,
				   cmplx* __restrict * __restrict padrevecs,
				   cmplx* __restrict * __restrict guessxt,
				   cmplx* __restrict * __restrict residxt,
					  temporalbasis* tbas );
  void cancelresiduals_generalized_eigenvector(int nbf,int ot, 
	cmplx* M1, cmplx* M2,cmplx* __restrict * __restrict guessxt,
	cmplx* __restrict * __restrict residxt,
	cmplx* __restrict * __restrict targetxt, temporalbasis* tbas); 
void cancelresiduals_generalized_eigenvector_Olap_preconditioner(int nbf,int ot, 
        cmplx* M1, cmplx* M2,cmplx* __restrict * __restrict guessxt,
	cmplx* __restrict * __restrict residxt,
	cmplx* __restrict * __restrict targetxt,temporalbasis* tbas);
  cmplx residualaction(int nbf,int ot,cmplx* __restrict * __restrict guessxt,
		       cmplx* __restrict * __restrict residxt);
  
  cmplx* __restrict 
    temporalvector_rhstarget(int ot, cmplx* __restrict levec,
			     cmplx* __restrict * __restrict psixt);
  cmplx* __restrict 
    temporalvectorprojection_revec(int ot, cmplx* __restrict levec,
			     cmplx* __restrict * __restrict psixt);
  cmplx* __restrict 
    temporalvectorprojection_levec(int ot, cmplx* __restrict revec,
			     cmplx* __restrict * __restrict psixt);
  void applycorrection_temporalevec(int ot,cmplx* __restrict dpsin, 
				    cmplx* __restrict vec,
				    cmplx* __restrict * __restrict guessxt,
				    cmplx* __restrict * __restrict residxt,
				    temporalbasis* tbas );
  void initialguess(int nmax, cmplx* __restrict psi0, 
		    cmplx* __restrict * __restrict psixt,
		    cmplx* __restrict * __restrict residxt,
		    temporalbasis* tbas);
  cmplx * __restrict * __restrict 
    deltawfxt_temporaleigenvectors(int nmax,int recursionlvl, 
				   int maxrecursionlvl, rl accuracygoal,
				   cmplx& stepaction,
				   cmplx* __restrict psi0, 
				   cmplx* __restrict * __restrict targetxt, 
				   temporalbasis* tbas);
  void timestep_temporaleigenvectors(rl tmin, rl tmax, potential* Pot,
				     pulse* Pls,wf* sourcewf, rl sourceEn,
				     rl accuracygoal, rl& newdt);

  //functions for use in fourier transform
  cmplx* fouriertransformgrid(rl xstart, rl xstop, int nftpts);
  cmplx* fouriertransformwf(rl xstart, rl xstop, int nftpts);
  cmplx* ftwf_shift(cmplx* __restrict wf_ft, int size, int nshift);
  cmplx* ftwf_recenter(cmplx* __restrict wf_ft, int size);
  cmplx* wf_fft_postprocess(rl xstart, rl xstop, int nftpts);
  rl* fftfrequencies(rl xstart, rl xstop, int nftpts);
  void printfft(str filename, rl xstart, rl xstop, int nftpts);
  void printfft(str filename, pulse* Pls,rl t,rl xstart, rl xstop, 
		int nftpts);
  
};



class wfxt{
 public:
  int nbf;
  int torder;
  globalbasis* gbas;
  cmplx* psixt;
  int* eltfirstindices;

  wfxt();
  wfxt(globalbasis* gbas_in);
  wfxt(wf* wf_in);
  ~wfxt();

  void propagate(pulse* Pls,rl tstart, rl tstop, rl dtstart,rl accuracygoal);
  void setupnextstep(pulse* Pls,rl t0, rl t1);

  void psi0setup(wf* psi0);
  void minimizeaction(cmplx* psi0);
  void minimizeaction();
  void wfxt_psi0xF0(cmplx* psi0);
  wf* psi_i();
  wf* psi_f();
  cmplx* psi_i_array();
  cmplx* psi_f_array();
  rl timestep(rl errortarget);

  
  
};
  
