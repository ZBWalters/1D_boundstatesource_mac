#include "./classes.h"
#include "./wf.h"


void wf::propagate(pulse* Pls,rl tstart, rl tstop, rl dtstart,rl accuracygoal){
  wfxt* tmpwfxt = new wfxt(this);
  tmpwfxt->propagate(Pls,tstart,tstop,dtstart,accuracygoal);
  delete [] psi;
  psi=tmpwfxt->psi_f_array();
  delete tmpwfxt;
}

//void wfxt::propagate(pulse* Pls,rl tstart,rl tstop,rl dtstart,rl accuracygoal){
//  //starting from a given wavefunction, take quality-controlled steps until
//  //final time has been reached.
//
//  // take first step
//  rl t0=tstart;
//  rl dt=dtstart;
//  rl t1=min(tstop,t0+dt);
//  gbas->tbasis->updateTmats(t0,t1,Pls);
//  minimizeaction();
//
//  // if first step hasn't reached the final time, enter a loop
//  while(t1<tstop){
//    t1=min(tstop,t0+dt);
//    setupnextstep(Pls,t0,t1);
//    minimizeaction();
//    rl dtnew=timestep(accuracygoal);
//    //dt=dtnew;
//    t0=t1;
//  }
//  
//}

void wfxt::propagate(pulse* Pls,rl tstart,rl tstop,rl dtstart,rl accuracygoal){
  //starting from a given wavefunction, take quality-controlled steps until
  //final time has been reached.

  cmplx* psii=psi_i_array();
  rl t0=tstart;
  rl t1=t0;
  rl dt=dtstart;

  // if first step hasn't reached the final time, enter a loop
  while(t1<tstop){
    t1=min(tstop,t0+dt);
    gbas->tbasis->updateTmats(t0,t1,Pls);
    wfxt_psi0xF0(psii);
    delete [] psii;
    minimizeaction();
    psii=psi_f_array();
    rl dtnew=timestep(accuracygoal);
    //dt=dtnew;
    t0=t1;
  }
  delete [] psii;
  
}

void wfxt::setupnextstep(pulse* Pls,rl t0, rl t1){
  gbas->tbasis->updateTmats(t0,t1,Pls);
  cmplx* psifarray=psi_f_array();
  cout << "psifarray\n";
  printmat(psifarray,nbf,1);
  wfxt_psi0xF0(psifarray);
  delete [] psifarray;
}
