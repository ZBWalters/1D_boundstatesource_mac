#include "./classes.h"
//#include "./potential.h"
//#include "./basis.h"
//#include "./legendrebasis.h"
//#include "./laguerrebasis.h"
#include "./derivedbasis.h"
#include "./mapintegrationtable.h"

rl** mapintegrationtable(legendrebasis* bas1, derivedbasis* bas2,
			 rl** integrationtable, int order){
  rl xmin=-1.;
  rl xmax=1.;
  rl ymin=bas2->xmin;
  rl ymax=bas2->xmax;

  rl dydx=(ymax-ymin)/(xmax-xmin);

  rl* intpts=new rl[order];
  rl* intwts=new rl[order];
  rl** rettable=new rl*[2];

  for(int i=0;i<order;i++){
    intpts[i]=xtoy(integrationtable[0][i],bas1,bas2);
    intwts[i]=integrationtable[1][i]*dydx;
  }
  rettable[0]=intpts;
  rettable[1]=intwts;
  return rettable;
  
}

rl xtoy(rl x,legendrebasis* bas1,derivedbasis* bas2){
  rl xmin=-1.;
  rl xmax=1.;
  rl ymin=bas2->xmin;
  rl ymax=bas2->xmax;

  return ymin+(x-xmin)/(xmax-xmin)*(ymax-ymin);
}
