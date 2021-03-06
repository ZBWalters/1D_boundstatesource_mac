#include "./classes.h"
//#include "./basis.h"
//#include "./legendrebasis.h"
//#include "./derivedbasis.h"
//#include "./globalbasis.h"
//#include "./pulse.h"
#include "./wf.h"

int main(){
  rl Pi=3.14159265358979323846;
  int ox=5;
  int ot=5;//15;//20;
  rl ecstheta=-0.7;//.7;//Pi/4.;//(Pi/4.);
  rl laguerrekappa=0.7;//ecstheta;//to get rid of single terms of kappa

  rl E0=.05;//0.001;
  rl w0=.0565;
  rl pulselength=1000.*Pi/w0;
  rl Qsource=1.;
  rl Qprop=1.;
  rl rho=1.;
  rl rhoprop=2.;

  pulse* Pls1=new pulse(E0,w0,pulselength);
  potential* potsource=new potential(Qsource,rho);
  potential* potprop=new potential(Qprop,rhoprop);

  //define tstart, tstop and interval between printed wfs
  rl dt=1.e-2;
  rl printstep=.01*Pi/w0;
  rl t0=-.5*Pi/w0;
  rl t1=t0+printstep;
  rl tstop=.5*Pi/w0;


//  int nelts=200;
//  int nbdyelts=20;
//  rl dx=20./nelts;
 

  rl Up=pow(E0,2.)/(4.*pow(w0,2.));
  rl En_kappa=max(1.*Up,4.);//max energy to use in determining element sizes
  cout << "Up, En_kappa\t"<<Up<<"\t"<<En_kappa<<"\n";
  rl deltaphi=Pi/5.;//Pi/5.;
  rl ecsbdy=10.;//100.;
  rl ecsbuffer=10.;//abs(log(1.e-6)/ecstheta);//100.;

  rl accuracygoal=1.e-6;//1.e-15;

  //cout << "setting up globalbasis\n";
  //globalbasis* gbas=new globalbasis(nelts,nbdyelts,dx,dt,ox,ot,laguerrekappa,ecstheta);
  //  gbas->ecssetup(nbdyelts,ecstheta);
  globalbasis* gbas_source=new globalbasis(potsource,Pls1,En_kappa,deltaphi, 
				    ecsbdy, ecsbuffer,dt,ox,ot, 
				    laguerrekappa,ecstheta);
  //use same inputs to ensure that the two globalbases are the same
  globalbasis* gbas_prop=new globalbasis(potsource,Pls1,En_kappa,deltaphi, 
				    ecsbdy, ecsbuffer,dt,ox,ot, 
				    laguerrekappa,ecstheta);
  
  gbas_source->setupVmats(potsource);
  gbas_source->updateHmats(potsource);
  gbas_prop->setupVmats(potprop);
  gbas_prop->updateHmats(potprop);
  

  //cmplx* psistart=initializearray(gbas->nbf);
  //cmplx (*startfun)(rl x)=exp(-pow(x,2));
  //cmplx* psistart=gbas->psi0setup(gauss);
  //cout << "psistart\n";
  //printmat(psistart,gbas->nbf,1);


  //cmplx* psinew=initializearray(gbas->nbf);

  wf* wf1=new wf(gbas_prop);
  wf* sourcewf=new wf(gbas_source);
  rl sourceEn=-.5;
 
//  //set up bound state as initial condition
  int nstate=0;//0;

  //wf1->boundstatesetup(nstate,sourceEn,t0);//set up bound state as initial condition
  //  wf1->psi_gaugechange(Pls1,t0,-1.);//change bound state to quasistatic state
  //use bound state as inhomogeneous source term
  sourcewf->boundstatesetup(nstate,sourceEn);//for use as source term
  
  ////interchange sourcewf->psi and wf1->psi
  cmplx* tmpwf=wf1->psi;
  wf1->psi=sourcewf->psi;
  //wf1->psi_gaugechange(Pls1,t0,-1.);//change bound state to quasistatic state
  sourcewf->psi=tmpwf;

 
  cout << "sourceEn\t"<<sourceEn<<"\n";

 
  //exchange source and propagated wfs
  //sourcewf->boundstatesetup(0,sourceEn);
  //cmplx* tmpwf=wf1->psi;
  //wf1->psi=sourcewf->psi;
  //sourcewf->psi=tmpwf;

//  //cout << "wf1\n";
//  //wf1->print();
//  ofstream wffile1 ("wf1.dat");
//  wffile1 << wf1->printstr();
//  wffile1.close();  newdt=dt*pow(abs(accuracygoal/error),(1./(ot-1.)));
//
//  
//  rl newdt=dt;
//  //wf1->timestep(0.,dt,pot1,Pls1,sourcewf,sourceEn,accuracygoal,newdt);
//  dt=1.;
//  wf1->propagate_iterative_inhomogenous(0.,1.,dt,Pls1,pot1,sourcewf,
//					sourceEn,accuracygoal);
//  ofstream wffile2 ("wf2.dat");
//  wffile2 << wf1->printstr();
//  wffile2.close();

  ofstream wffile1("wf0.dat");
  wffile1 << wf1->printstr();
  wffile1.close();
  wffile1.clear();

  ofstream file2("norm_vs_time.dat");

  rl newdt=dt;
  //wf1->timestep(0.,dt,pot1,Pls1,sourcewf,sourceEn,accuracygoal,newdt);
  int printindx=0;
  while(t0<tstop){
    //mtrace();
    wf1->propagate_iterative_inhomogenous(t0,t1,dt,Pls1,potprop,sourcewf,
					  sourceEn,accuracygoal);
    //muntrace();
//    wf1->propagate_iterative(Pls1,pot1,ot,t0,t1,dt,accuracygoal);
    printindx++;
    wffile1.open(("wf"+to_string(printindx)+".dat").c_str());
    wffile1 << wf1->printstr();
    wffile1.close();
    wffile1.clear();
    file2<<t1<<"\t"<<wf1->psinorm()<<"\n";
    //wf1->printfft("fft"+to_string(printindx)+".dat",Pls1,t1,-ecsbdy*.99,ecsbdy*.99,4096);
    //wf1->printfft("fft"+to_string(printindx)+".dat",-ecsbdy*.99,ecsbdy*.99,4096);
    t0=t1;
    t1=t0+printstep;
  }
  file2.close();


  //cout << "psif\n";
  //printmat(wf1->psi,wf1->nbf,1);


  //old way
    //Propagate wf1 for a given time & plot result
//  wf1->propagate(Pls1,0.,1.,.005,1.e-6);
//  wf* wf2=wfxt1->psi_f();
//  ofstream wffile2 ("wf2.dat");
//  wffile2 << wf1->printstr();
//  wffile2.close();



  delete gbas_source;
  delete gbas_prop;
  delete potsource;
  delete potprop;
  delete sourcewf;
  delete wf1;
  //  delete wf2;
  //delete wfxt1;
  delete Pls1;

}
