#include "classes.h"

void printmat(rl* mat, int n1, int n2){
  //Print out matrix
  for(int i=0;i<n1;i++){
    cout << "i "<< i << "\t";
    for(int j=0;j<n2;j++){
      cout << mat[arrayindx(i,j,n1,n2)] << "\t";
    }
    cout << "\n";
  }
  cout << "\n";
}

void printmat(cmplx* mat, int n1, int n2){
  //Print out matrix
  for(int i=0;i<n1;i++){
    cout << "i "<< i << "\t";
    for(int j=0;j<n2;j++){
      cout << mat[arrayindx(i,j,n1,n2)] << "\t";
    }
    cout << "\n";
  }
  cout << "\n";
}

void printmat_scipy(cmplx* mat, int n1, int n2){
  //Print out matrix
  cout << "[";
  for(int i=0;i<n1;i++){
    cout << "[";
    for(int j=0;j<n2-1;j++){
      cout << real(mat[arrayindx(i,j,n1,n2)])<<"+1j*" << imag(mat[arrayindx(i,j,n1,n2)])<<",\t";
    }
    cout << real(mat[arrayindx(i,n2-1,n1,n2)])<<"+1j*" << imag(mat[arrayindx(i,n2-1,n1,n2)])<<"],";
    cout << "\n";
  }
  cout << "]";
  cout << "\n";
}

void printmat_scipy(rl* mat, int n1, int n2){
  //Print out matrix
  cout << "[";
  for(int i=0;i<n1;i++){
    cout << "[";
    for(int j=0;j<n2-1;j++){
      cout << mat[arrayindx(i,j,n1,n2)]<<",\t";
    }
    cout << mat[arrayindx(i,n2-1,n1,n2)]<<"],";
    cout << "\n";
  }
  cout << "]";
  cout << "\n";
}

void printmat_banded(cmplx* mat,int n1, int kl, int ku, int n2){
  //print out the nonzero band of a matrix using lapack's banded matrix
  //format, as given in the man page for zgbsv
  for(int j=0;j<n2;j++){
    int imin=max(0,j-ku);
    int imax=min(n1-1,j+kl);
    cout<< imin<< " "<<j<<" [";
    for(int i=imin;i<=imax;i++){
      cout<< mat[bandarrayindx(i,j,ku,kl,n2)]<< "\t";
    }
    cout << "]\n";
  }
}

void printmat_formatted(str filename, cmplx* mat, int n1, int n2){
  //print out matrix in form suitable for plotting
  const char* cfilename=filename.c_str();
  ofstream file1 (cfilename);
  file1<<matstr_formatted(mat,n1,n2);
  file1.close();
}

str matstr_formatted(cmplx* mat, int n1,int n2){
  //Print out matrix
  str retstr="";
  for(int i=0;i<n1;i++){
    //retstr+= "i "+to_string(i) + "\t";
    retstr+= to_string(i) + "\t";//i index not useful for plotting
    for(int j=0;j<n2;j++){
      retstr+=to_string(real(mat[arrayindx(i,j,n1,n2)]))+"\t"+to_string(imag(mat[arrayindx(i,j,n1,n2)])) + "\t";
    }
    retstr+= "\n";
  }
  retstr+= "\n";
  return retstr;
}

str matpairstr_formatted(rl* xmat, cmplx* zmat, int n){\
  str matstr="";
  for(int i=0;i<n;i++){
    matstr+=to_string(xmat[i])+"\t"+to_string(real(zmat[i]))+"\t"+
      to_string(imag(zmat[i]))+"\n";
  }
  return matstr;
}

void printmatpair(str filename,rl* xmat, cmplx* zmat, int n){
  const char* cfilename=filename.c_str();
  ofstream file1 (cfilename);
  file1<<matpairstr_formatted(xmat,zmat,n);
  file1.close();
}
