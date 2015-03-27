from numpy import *
from scipy import *
from scipy.integrate import quad


wL=.0565#laser frequency
F=0.05#field strength
Upon=(F**2)/(4*wL**2)#ponderomotive potential
Ip=.669#0.5#Ionization potential
gammaKeldysh=(Ip/(2.*Upon))**.5#Keldysh parameter
tau=gammaKeldysh/wL#tunneling time
V0=abs(F/wL)

def keldyshparameter(F,Ip,wL):
    Up=(F**2.)/(4.*wL**2.)
    gammaK=(abs(Ip)/(2.*Up))**.5
    return gammaK
 

def Tauion(gkel):
    return arcsinh(gkel)/wL
def Tion(gkel):
    return 1j*Tauion(gkel)
def RTion(gkel):
    return real(Tion(gkel))


def amp1(Ip,gkel):
    return exp(1j*Ip*Tion(gkel))

def V(pconjg,t):
    return pconjg-(F/wL)*sin(wL*t)

def KE_r(tr,pconj):
    return 0.5*V(pconjg,tr)**2.

def KE_im(ti,pconjg):
    return 0.5*V(pconjg,ti*1j)**2.
def rKE_im(ti,pconjg):
    return real(KE_im(ti,pconjg))
def imKE_im(ti,pconjg):
    return imag(KE_im(ti,pconjg))


def amp2(pconjg,gkel):
    t1=Tauion(gkel)
    rretval1=quad(rKE_im,t1,0,args=(pconjg))[0]
    imretval1=quad(imKE_im,t1,0,args=(pconjg))[0]
    phase=-1j*(rretval1+imretval1*1j)
    return exp(1j*phase)

def amp(pconjg,gkel):
    return amp1(Ip,gkel)*amp2(pconjg,gkel)

##################################
###integrate v(t) to find positions at time of measurement
def vel_tim(t_im,pconjg,gkel):
    return 1j*(pconjg-F/wL*sin(wL*(RTion(gkel)+1j*t_im)))
def rvel_tim(t_im,pconjg,gkel):
    return real(vel_tim(t_im,pconjg,gkel))
def imvel_tim(t_im,pconjg,gkel):
    return imag(vel_tim(t_im,pconjg,gkel))


def vel_tr(tr,pconjg,gkel):
    return pconjg-F/wL*sin(wL*(RTion(gkel)+tr))
def rvel_tr(tr,pconjg,gkel):
    return real(vel_tr(tr,pconjg,gkel))
def imvel_tr(tr,pconjg,gkel):
    return imag(vel_tr(tr,pconjg,gkel))

def Xion(pconjg,gkel):
    tion=Tauion(gkel)
    retval1=quad(rvel_tim,tion,0,args=(pconjg,gkel))[0]
    retval2=quad(imvel_tim,tion,0,args=(pconjg,gkel))[0]
    return retval1+1j*retval2

def Xm(pconjg,gkel,tm):
    retval1=quad(rvel_tr,RTion(gkel),tm,args=(pconjg,gkel))[0]
    retval2=quad(imvel_tr,RTion(gkel),tm,args=(pconjg,gkel))[0]
    return retval1+1j*retval2

def X(pconjg,gkel,tm):
    return real(Xm(pconjg,gkel,tm))#real(Xion(pconjg,gkel)+Xm(pconjg,gkel,tm))

########################################################

deltap=.01*V0
ptable=arange(-V0,V0,deltap)

xtable=ptable*0.
amp1table=ptable*0j
amp2table=ptable*0j
amptable=ptable*0j

tplot=0.5*pi/wL

for i in range(len(ptable)):
    xtable[i]=X(ptable[i],gammaKeldysh,tplot)
    amp1table[i]=amp1(Ip,gammaKeldysh)
    amp2table[i]=amp2(ptable[i],gammaKeldysh)
    amptable[i]=amp(ptable[i],gammaKeldysh)

f=open('PPT.dat','w')
for i in range(len(xtable)):
    f.write(str(real(xtable[i]))+"\t"+str(real(amptable[i]))+"\t"+str(imag(amptable[i]))+"\n")
f.close()
