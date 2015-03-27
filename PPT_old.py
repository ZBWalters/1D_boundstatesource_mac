from numpy import *
from scipy import *
from scipy.integrate import quad


wL=.0565#laser frequency
F=-0.05#field strength
Upon=(F**2)/(4*wL**2)#ponderomotive potential
Ip=0.5#Ionization potential
gammaKeldysh=(Ip/(2.*Upon))**.5#Keldysh parameter
tau=gammaKeldysh/wL#tunneling time
V0=abs(F/wL)

def Tauion(pconjg,gkel):
    return arcsinh(gkel-1j*pconjg/V0)/wL

def Tion(pconjg,gkel):
    return Tauion(pconjg,gkel)*1j

def RTion(pconjg,gkel):
    return real(Tion(pconjg,gkel))

def Tbirth(pconjg):
    return arcsin(pconjg/V0)/wL


def S01(pconjg,gkel,tm):
    retval1=Ip*(Tion(pconjg,gkel))
    return retval1

def S02(pconjg,gkel,tm):
    rretval2=quad(rfunc2,Tauion(pconjg,gkel),0,args=(pconjg,gkel))[0]
    imretval2=quad(imfunc2,Tauion(pconjg,gkel),0,args=(pconjg,gkel))[0]
    retval2=-1j*(rretval2+1j*imretval2)
    return retval2

def S03(pconjg,gkel,tm):
    rretval3=quad(rfunc3,RTion(pconjg,gkel),tm,args=(pconjg,gkel))[0]
    imretval3=quad(imfunc3,RTion(pconjg,gkel),tm,args=(pconjg,gkel))[0]
    retval3=rretval3+1j*imretval3
    return retval3

def S0(pconjg,gkel,tm):
    return S01(pconjg,gkel,tm)+S02(pconjg,gkel,tm)+S03(pconjg,gkel,tm)

def func2(tau,pconjg,gkel):
    return (0.5+0j)*(pconjg-(F/wL)*
                          sin(wL*(RTion(pconjg,gkel)+1j*tau)))**2.
def rfunc2(tau,pconjg,gkel):
    return real(func2(tau,pconjg,gkel))
def imfunc2(tau,pconjg,gkel):
    return imag(func2(tau,pconjg,gkel))
    
def func3(t,pconjg,gkel):
    return (0.5+0j)*(pconjg-(F/wL)*sin(wL*t))**2.
def rfunc3(t,pconjg,gkel):
    return real(func3(t,pconjg,gkel))
def imfunc3(t,pconjg,gkel):
    return imag(func3(t,pconjg,gkel))

def amplitude(pconjg,gkel,tm):
    return exp(1j*S0(pconjg,gkel,tm))


def expfactor(Ip,F):
    return exp((-2.*(2.*Ip)**(1.5))/(3.*abs(F)))

def B(pconjg,gkel):
    tion=Tion(pconjg,gkel)
    eps=eps_paper(pconjg,gkel,tion)
    Ftmp=F*sin(wL*tion)
    return B_paper(Ftmp,Ip,eps)

def B_paper(F0,Ip,eps):#from Eq 14 in Popruzhenko paper
    return 1j/sqrt(2.*pi*F0)*exp((-2.*sqrt(2)/3.)*((Ip+eps)**1.5)/F0)

def eps_paper(pconjg,gkel,t):
    return .5*(pconjg-V0*sin(wL*t))

###integrate v(t) to find positions at time of measurement
def vel_tim(t_im,pconjg,gkel):
    return 1j*(pconjg-F/wL*sin(wL*(RTion(pconjg,gkel)+1j*t_im)))
def rvel_tim(t_im,pconjg,gkel):
    return real(vel_tim(t_im,pconjg,gkel))
def imvel_tim(t_im,pconjg,gkel):
    return imag(vel_tim(t_im,pconjg,gkel))


def vel_tr(tr,pconjg,gkel):
    return pconjg-F/wL*sin(wL*(RTion(pconjg,gkel)+tr))
def rvel_tr(tr,pconjg,gkel):
    return real(vel_tr(tr,pconjg,gkel))
def imvel_tr(tr,pconjg,gkel):
    return imag(vel_tr(tr,pconjg,gkel))

def Xion(pconjg,gkel):
    retval1=quad(rvel_tim,Tauion(pconjg,gkel),0,args=(pconjg,gkel))[0]
    retval2=quad(imvel_tim,Tauion(pconjg,gkel),0,args=(pconjg,gkel))[0]
    return retval1+1j*retval2

def Xm(pconjg,gkel,tm):
    retval1=quad(rvel_tr,RTion(pconjg,gkel),tm,args=(pconjg,gkel))[0]
    retval2=quad(imvel_tr,RTion(pconjg,gkel),tm,args=(pconjg,gkel))[0]
    return retval1+1j*retval2

def X(pconjg,gkel,tm):
    return real(Xion(pconjg,gkel)+Xm(pconjg,gkel,tm))

######make tables & plot
deltap=.01
pconjgtable=arange(-V0,V0,deltap)
amptable=pconjgtable*0j
xtable=pconjgtable*0j
Stable=pconjgtable*0j
Stable1=pconjgtable*0j
Stable2=pconjgtable*0j
Stable3=pconjgtable*0j
Ttable=pconjgtable*0j
tplot=.5*pi/wL
for i in range(len(pconjgtable)):
    pconjg=pconjgtable[i]
    xtable[i]=X(pconjg,gammaKeldysh,tplot)
    amptable[i]=abs(amplitude(pconjgtable[i],gammaKeldysh,tplot))
    Stable[i]=S0(pconjgtable[i],gammaKeldysh,tplot)
    Stable1[i]=S01(pconjgtable[i],gammaKeldysh,tplot)
    Stable2[i]=S02(pconjgtable[i],gammaKeldysh,tplot)
    Stable3[i]=S03(pconjgtable[i],gammaKeldysh,tplot)
    Ttable[i]=Tion(pconjgtable[i],gammaKeldysh)
