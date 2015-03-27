#uses formulas from Attosecond correlation dynamics during
#electron tunnelling from molecules, Z. Walters & O. Smirnova 2010
from numpy import *
from scipy import *
from scipy.integrate import quad

wL=.0565#laser frequency
F=0.05#field strength
Upon=(F**2)/(4*wL**2)#ponderomotive potential
Ip=.669#0.5#Ionization potential
gammaKeldysh=(Ip/(2.*Upon))**.5#Keldysh parameter
tau=gammaKeldysh/wL#tunneling time
V0=-F/wL#abs(F/wL)


def Tion(kf,Ip):
    ionizationphase=arcsin((1j*sqrt(2.*Ip)-kf)/abs(V0))
    ionizationtime=ionizationphase/wL
    return ionizationtime

def vecA(t):
    return V0*sin(wL*t+0j)

def vel(t,kf):
    return kf+vecA(t)

def S(t,kf):
    return .5*vel(t,kf)**2.

def rS(t,kf):
    return S(t,kf).real

def imS(t,kf):
    return S(t,kf).imag

def rS_barrier(ti,tr,kf):
    return rS(tr+1j*ti,kf)

def imS_barrier(ti,tr,kf):
    return imS(tr+1j*ti,kf)

def rvel(t,kf):
    return vel(t,kf).real

def imvel(t,kf):
    return vel(t,kf).imag

def rvel_barrier(tim,tr,kf):#velocity while under the barrier
    return (vel(tr+1j*tim,kf)).real

def imvel_barrier(tim,tr,kf):#velocity while under the barrier
    return (vel(tr+1j*tim,kf)).imag

def rX(t,kf,Ip):#real part of x position
    t0=Tion(kf,Ip)
    rt0=t0.real
    it0=t0.imag
    
    #real axis just below time of ionization
    t1=rt0

    #movement under the barrier
    #real part of delta x = 1j*imvel
    x_tunnel=-quad(imvel_barrier,it0,0,args=(rt0,kf))[0]

    x_continuum=quad(rvel,rt0,real(t),args=(kf))[0]
    return x_tunnel+x_continuum

def tunnelphase(kf,Ip):
    t0=Tion(kf,Ip)
    rt0=t0.real
    it0=t0.imag

    rtunnelphase=quad(rS_barrier,it0,0,args=(rt0,kf))[0]
    imtunnelphase=quad(imS_barrier,it0,0,args=(rt0,kf))[0]

    tunnelphase=rtunnelphase+1j*imtunnelphase
    tunnelphase=tunnelphase*-1j
    return tunnelphase

def continuumphase(t,kf,Ip):
    t0=Tion(kf,Ip)
    rt0=t0.real

    rcontinuumphase=quad(rS,rt0,t,args=(kf))[0]
    imcontinuumphase=quad(imS,rt0,t,args=(kf))[0]

    return rcontinuumphase+1j*imcontinuumphase

def boundstatephase(kf,Ip):
    return Ip*Tion(kf,Ip)

def totalphase(t,kf,Ip):
    return tunnelphase(kf,Ip)+continuumphase(t,kf,Ip)+boundstatephase(kf,Ip)

def amplitude(t,kf,Ip):
    return exp(1j*totalphase(t,kf,Ip))


######################################

tplot=.5*pi/wL

deltak=.01*V0
ktable=arange(-V0,V0,deltak)

xtable=ktable*0.
amptable=ktable*0j

for i in range(len(ktable)):
    xtable[i]=rX(tplot,ktable[i],Ip)
    amptable[i]=amplitude(tplot,ktable[i],Ip)

f=open('PPT.dat','w')
for i in range(len(xtable)):
    f.write(str(xtable[i].real)+"\t"+str(amptable[i].real)+"\t"+str(amptable[i].imag)+"\n")
f.close()
