#include "./classes.h"
#include "./pulse.h"
//typedef double rl;

class pulse;

pulse::pulse(){
  E0=.0;//defaultE0;
  omega=.0565;//defaultomega;
  pulselength=1000.;//defaultpulselength;
  A0=-E0/omega;
}

pulse::pulse(rl in_E0, rl in_omega){
  E0=in_E0;
  omega=in_omega;
  pulselength=1000;//defaultpulelength;
  A0=-E0/omega;
}

pulse::pulse(rl in_E0, rl in_omega, rl in_pulselength){
  E0=in_E0;
  omega=in_omega;
  pulselength=in_pulselength;
  A0=-E0/omega;
}




rl pulse::Ez(rl t){
  rl retval=E0*cos(omega*t);//*exp(-pow(t/pulselength,2));
//  rl retval=E0*cos(omega*t)*exp(-pow(t/pulselength,2))-
//    (2.*E0*t/(pow(pulselength,2.)*omega))*sin(omega*t)*exp(-pow(t/pulselength,2));
  return retval;
}

rl pulse::Az(rl t){
  //cout << "Az t "<< t <<"\n";
  rl retval=A0*sin(omega*t);
  //rl retval= A0*sin(omega*t)*exp(-pow(t/pulselength,2));
  return retval;
}
