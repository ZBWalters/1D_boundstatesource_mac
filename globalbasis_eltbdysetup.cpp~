#include "./classes.h"
#include "./globalbasis.h"

//methods in this file use wkb to find element sizes such that a wavefunction
//of energy En accumulates a phase of approximately deltaphi within a single
//element.  Starting at x=0, these elements will be joined together until the
//desired xmax is reached.  In order to conform to the needs of exterior
//complex scaling, one element boundary will be located at ecsbdy, while the
//final, infinite range element will be required to start at xmin>
//ecsbdy+ecsbuffer, to give the travelling waves some time to decay before
//reaching the infinite range element.

rl globalbasis::nextdx(potential* pot,rl x, rl En, rl deltaphi){
  //find size of basis element such that sqrt(2(E-V))*dx=deltaphi

  rl dx=0.0;
//  if(x==0.){
//    //dx=.01;
//    dx=
//  }
//  else{
    //rl kappa=sqrt(2.*abs(En+abs(pot->V(x))));
  rl V=max(abs(pot->V(x)),abs(1./x));
  rl kappa=sqrt(2.*(abs(En)+V));
  dx=deltaphi/kappa;
    //  }
  //cout << "next dx "<<dx<<"\n";
  return dx;
}

rl globalbasis::nextdx(potential* pot,pulse* Pls,rl x,rl ecsbdy, rl ecstheta,rl En, rl deltaphi){
  //find size of basis element such that sqrt(2(E-V))*dx=deltaphi

  rl dx=0.;
//  if(x==0.){
//    dx=.01;
//    //dx=pow((deltaphi/(2.*sqrt(2.))),2);
//    //dx=min(dx,pot->rc);
//  }
//  else{
    //rl kappa=sqrt(2.*abs(En+abs(pot->V(x))));
    rl E0=Pls->E0;
    rl V=abs(pot->V(x))+abs(E0*x);//abs(pot->V(x));//abs(pot->V(x))+abs(E0*x);//max(abs(pot->V(x))+abs(E0*x),abs(1./x));
    rl kappa=sqrt(2.*(abs(En)+V));
    dx=deltaphi/kappa;
//  }

    if(abs(x)>abs(ecsbdy)){
      dx*=cos(2.*ecstheta);//.1;//sqrt(2.)/2.;//sin(Pi/4.);
  }
  //cout << "next dx "<<dx<<"\n";
  return dx;
}


rl globalbasis::nextdx(potential* pot,rl x,rl ecsbdy, rl ecstheta,rl En, rl deltaphi){
  //find size of basis element such that sqrt(2(E-V))*dx=deltaphi

  rl dx=0.;
  //if(x==0.){
  //  dx=.1;
  //  //dx=pow((deltaphi/(2.*sqrt(2.))),2);
  //  //dx=min(dx,pot->rc);
  //}
  //else{
    //rl kappa=sqrt(2.*abs(En+abs(pot->V(x))));
    rl V=abs(pot->V(x));//max(abs(pot->V(x)),abs(1./x));
    rl kappa=sqrt(2.*(abs(En)+V));
    dx=deltaphi/kappa;
    //}

    if(abs(x)>abs(ecsbdy)){
      dx*=cos(2.*ecstheta);//.1;//sqrt(2.)/2.;//sin(Pi/4.);
  }
  //cout << "next dx "<<dx<<"\n";
  return dx;
}




rl* globalbasis::elementboundaries_zerobased(potential* pot,pulse* Pls, rl En, 
					     rl deltaphi, rl ecsbdy, 
					     rl ecsbuffer, rl ecstheta, 
					     int &nbdy, 
					     int &necselts, 
					     bool returnarray){
  //if retarray= true, return a non-null array of boundary points.  Otherwise,
  //return null (since nbdy is passed by reference, this will simply calculate
  //the size of the necessary array.
  
  rl ecsbdy_out=ecsbdy+ecsbuffer;

  rl* retarray=0;
  if(returnarray){
    retarray=new rl[nbdy];
  }
  nbdy=1;
  necselts=0;
  rl x1=0.;
  rl x2=0.;
  if(returnarray){
    retarray[nbdy-1]=x1;
  }
  
  while(x1<ecsbdy_out){
    rl dx=nextdx(pot,Pls,x1,ecsbdy,ecstheta,En,deltaphi);
    x2=x1+dx;
    if((x1<ecsbdy) and (x2>ecsbdy)){
      x2=ecsbdy;
    }
    nbdy++;
    if(x1>=ecsbdy){
      necselts++;
    }
    x1=x2;
    if(returnarray){
      retarray[nbdy-1]=x1;
    }
  }
  //necselts++;
  //nbdy++;//size of retarray is equal to index of last point + 1;
  return retarray;
}

rl* globalbasis::elementboundaries_zerobased(potential* pot, pulse* Pls,rl En, 
					     rl deltaphi, rl ecsbdy, 
					     rl ecsbuffer, rl ecstheta, 
					     int &nbdy,
					     int &necselts){
  //calls elementboundaries_zerobased once with returnarray set to false (to
  //calculate correct array size), then once with returnarray set to true, to
  //return the correctly sized array;
  bool returnarray=false;
  rl* retval=elementboundaries_zerobased(pot,Pls,En,deltaphi,ecsbdy,ecsbuffer,
					 ecstheta,nbdy,necselts,returnarray);
  returnarray=true;
  //cout <<"nbdy after 1st call "<<nbdy<<"\n";
  //cout <<"necselts after 1st call "<<necselts<<"\n";
  retval=elementboundaries_zerobased(pot,Pls,En,deltaphi,ecsbdy,ecsbuffer,
				     ecstheta,nbdy,necselts,returnarray);
  return retval;
}
