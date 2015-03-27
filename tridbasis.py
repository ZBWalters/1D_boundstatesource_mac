from numpy import *
from scipy import *
from scipy.linalg import eig

def triD_U(n,m):
    retval=0.
    nmin=min(n,m)
    nmax=max(n,m)
    deltan=nmax-nmin
    if(nmin==0):
        if(nmax==0):
            retval=2./3.
        if(nmax==1):
            retval=1./3.
        if(nmax==2):
            retval=-1./6.
    else:
        if(deltan==0):
            retval=(2.*nmax)/(4.*(nmax**2)-1.)
        if(deltan==1):
            retval=(2./(2.*nmin+1.))/4.
    return retval

def triD_Q(n,m):
    retval=0.
    if(n>m):
        retval=1.
    if(n==m):
        retval=1./2.
    if(m==0):
        if(n==1):
            retval=1./2.
        else:
            retval=0.
    if(n==0):
        if(m==0):
            retval=-1./2.
        if(m==1):
            retval=-1./2.
    return retval

def triD_Qinv(n):
    Q=triD_Qmat(n)[1:,1:]
    return linalg.inv(Q)

def triD_QinvU(n):
    return dot(triD_Qinv(n),triD_Umat(n)[1:,1:])

def triD_UQinv(n):
    return dot(triD_Umat(n)[1:,1:],triD_Qinv(n))

def triD_UQinv_uppertriangle(n):
    retmat=triD_UQinv(n)
    for i in range(n-1):
        for j in range(i+1,n-1):
            retmat[j,i]=0.
    return retmat


def triD_QinvU_uppertriangle(n):
    retmat=triD_QinvU(n)
    for i in range(n-1):
        for j in range(i+1,n-1):
            retmat[j,i]=0.
    return retmat



def triD_QinvU_lowertriangle(n):
    retmat=triD_QinvU(n)
    for i in range(n-1):
        for j in range(i+1,n-1):
            retmat[i,j]=0.
    return retmat


def utsolveeigs(n):
    mat=eye(n-1)-dot(triD_QinvU(n),linalg.inv(triD_QinvU_uppertriangle(n)))
    return linalg.eig(mat)[0]

def triD_Umat(n):
    retarray=zeros((n,n))
    for i in range(n):
        for j in range(n):
            retarray[i,j]=triD_U(i,j)
    return retarray

def triD_Umat_legendre(n):
    retarray=zeros((n+1,n+1))
    for i in range(n):
        for j in range(n):
            retarray[i,j]=triD_U(i,j)
    retarray[0,n]=1.
    retarray[n,0]=1.
    return retarray

def triD_legendre_vec(n):
    lhsmat=triD_Umat_legendre(n)
    rhsvec=zeros(n+1)
    rhsvec[0]=1.
    return linalg.solve(lhsmat,rhsvec)

def triD_Umat_lstsq(n):
    retmat=triD_Qmat(n)
    return retmat[:,1:]

def triD_lstsq_vec(n):
    lhsmat=triD_Umat_lstsq(n)
    rhsvec=zeros(n)
    rhsvec[0]=1.
    lhsmat=dot(triD_sqrtU(n),lhsmat)
    rhsvec=dot(triD_sqrtU(n),rhsvec)
    return linalg.lstsq(lhsmat,rhsvec)

def triD_lstsq_Udotvec(n):
    smallvec=triD_lstsq_vec(n)[0]
    bigvec=append([0.],smallvec)
    return dot(triD_Umat(n),bigvec)

def triD_Uinv(n):
    return linalg.inv(triD_Umat(n))


def triD_Qmat(n):
    retarray=zeros((n,n))
    for i in range(n):
        for j in range(n):
            retarray[i,j]=triD_Q(i,j)
    return transpose(retarray)

def triD_IC(n):
    icvec=zeros(n)
    icvec[0]=1.
    icvec[1]=0.
    #print "icvec",icvec
    residvec=dot(triD_Umat(n),icvec)
    #print "residvec",residvec
    for i in range(1,n):
        tmpvec=zeros(n)
        tmpvec[i]=residvec[i]
        deltaicvec=dot(triD_Uinv(n),tmpvec)
        for j in range(n):
            icvec[j]=icvec[j]-deltaicvec[j]
    tmpval=icvec[0]
    icvec=icvec/tmpval
    return icvec

def vec(n,m):
    retvec=zeros(n)
    retvec[m]=1.
    return retvec

def QUeig(ot):
    QR=triD_Qmat(ot)[1:,1:]
    UR=triD_Umat(ot)[1:,1:]
    return eig(QR,UR)

def triD_basismat(n):
    basismat=zeros((n,n))
    basismat[0,0]=1.
    basismat[0,1]=-1.
    for m in range(1,n):
        basismat[m,m]=1.
        basismat[m,m-1]=1.
    basismat=basismat/2.
    return basismat


def triD_sqrtU(n):
    legendreweights=zeros((n,n))
    for m in range(0,n):
        legendreweights[m,m]=sqrt(2./(2.*m+1.))
    basismat=triD_basismat(n)
    #return dot(transpose(basismat),dot(legendreweights,basismat))
    return dot(basismat,dot(legendreweights,transpose(basismat)))

def triD_ic0_vec(n):
    Utmp=triD_Umat(n)
    vec=zeros(n)
    vec[0]=1.
    retvec=linalg.solve(Utmp,vec)
    norm=dot(transpose(retvec),dot(Utmp,retvec))
    return retvec/sqrt(norm)

