//apply corrections in the same way that the linear system is set up
void wf::applydeltapsi_Aresidual(int ncorr,int nmax, cmplx* __restrict deltapsi,
				 cmplx Aratio, cmplx Asqratio, cmplx* __restrict * __restrict residxt, 
		       temporalbasis* tbas ){
  cmplx II=cmplx(0.,1.);
  rl ecstheta=gbas->ecstheta;
  int ot=tbas->order;
  //apply deltapsi vector to entry n of guessxt, adjust residuals by
  //L.deltapsi

  cout << "inside applydeltapsi\n";
  cout <<"umat\n";
  printmat(tbas->umat,ot,ot);
  cout <<"qmat\n";
  printmat(tbas->qmat,ot,ot);

  //apply deltapsi vector to entry n of guessxt
  //for(int i=0;i<nbf;i++){
  //  guessxt[ncorr][i]+=deltapsi[i];
  //}
  int inc=1;
  cmplx one=1.;

  for(int eltnum=0;eltnum<gbas->nelts;eltnum++){
    derivedbasis* tmpbasis=gbas->elementbases[eltnum];
    int tmporder=tmpbasis->order;
    int eltindx=gbas->eltfirstindices[eltnum];
    cmplx* __restrict tmpo=tmpbasis->omat;
    cmplx* __restrict tmpnabla=tmpbasis->nablamat;
    cmplx* __restrict tmptransnabla=transpose(tmpnabla,tmporder,tmporder);


    //cout << "tmpnabla\n";
    //printmat(tmpnabla,tmporder,tmporder);
    //cout << "tmpnablasq\n";
    //printmat(tmpnablasq,tmporder,tmporder);


    cmplx alpha=1.;
    cmplx beta=0.;

    //find product of Hamiltonian & overlap arrays within temporary basis
    cmplx* __restrict tmpodeltapsi=tmpbasis->Xmatpsi(tmpbasis->omat,
						&deltapsi[eltindx]);
        cmplx* __restrict tmpnabladeltapsi=tmpbasis->Xmatpsi(tmpbasis->nablamat,
						&deltapsi[eltindx]);

    cmplx* __restrict tmptransposenabladeltapsi=tmpbasis->
      Xmatpsi(tmptransnabla,&deltapsi[eltindx]);
    
    //sum L.psi for current element & apply to residxt
    cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)]-
      Aratio*tbas->umat[arrayindx(n,ncorr,ot,ot)];
      if(tmpbasis->ecsflag){
	coeff=0.5*II*Aval*exp(II*ecstheta);
      }
      else{
	coeff=0.5*II*Aval;
      }
      cblas_zaxpy(tmporder,&coeff,tmpnabladeltapsi,inc,tmpdeltaresid,inc);

      if(tmpbasis->ecsflag){
	coeff=-0.5*II*Aval*exp(-II*ecstheta);
      }
      else{
	coeff=-0.5*II*Aval;
      }
      cblas_zaxpy(tmporder,&coeff,tmptransposenabladeltapsi,inc,
		  tmpdeltaresid,inc);
      cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)]-
	Asqratio*tbas->umat[arrayindx(n,ncorr,ot,ot)];
      if(tmpbasis->ecsflag){
	coeff=0.5*Asqval*exp(2.*II*ecstheta);
      }
      else{
	coeff=0.5*Asqval;
      }
      cblas_zaxpy(tmporder,&coeff,tmpodeltapsi,inc,
		  tmpdeltaresid,inc);

      //apply tmpdeltaresid to residxt[n]
      cblas_zaxpy(tmporder,&one,tmpdeltaresid,inc,&residxt[n][eltindx],inc);

      delete [] tmpdeltaresid;
    }

    
    delete [] tmpodeltapsi;
    delete [] tmpnabladeltapsi;
    delete [] tmptransnabla;
    delete [] tmptransposenabladeltapsi;

    

  }
}










cmplx Aval=tbas->Amat[arrayindx(n,ncorr,ot,ot)];
      if(tmpbasis->ecsflag){
	coeff=0.5*II*Aval*exp(II*ecstheta);
      }
      else{
	coeff=0.5*II*Aval;
      }


if(tmpbasis->ecsflag){
	coeff=-0.5*II*Aval*exp(-II*ecstheta);
      }
      else{
	coeff=-0.5*II*Aval;
      }

cmplx Asqval=tbas->Asqmat[arrayindx(n,ncorr,ot,ot)];
      if(tmpbasis->ecsflag){
	coeff=0.5*Asqval*exp(2.*II*ecstheta);
      }
      else{
	coeff=0.5*Asqval;
      }
