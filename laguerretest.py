#simple script to work out how a function should be expanded in terms of
#laguerre polynomials.  the tricky detail is that the gauss laguerre weights
#already include the factor of exp(-x) needed to give orthogonality, so that
#\sum_{i} L_n(x_i) L_m(x_i) w_i = delta_{n,m} In the code, I would like to use
#basis functions consisting of Laguerre polynomials times exponentials, or
#bf_n(x)=L_n(x) exp(-x/2), so that \int_{0}^{\infty} bf_n(x)
#bf_m(x)=delta_{n,m}.  With respect to this basis, expansion coefficients are
#found via
# A_n=\sum_i f(x_i) exp(x_i/2) bf_n(x_i) w_i


import scipy
import numpy
from pylab import *

from numpy import exp,zeros,arange
from scipy.special import laguerre,airy
from scipy.special.orthogonal import la_roots

def twofunctionintegral(n,m):
    order=2*max(n,m)
    [glpts,glwts]=la_roots(order,0)
    print "glpts ",glpts
    print "glwts ",glwts
    sum=0.;
    for i in range(len(glpts)):
        x=glpts[i]
        Ln=laguerre(n)(x)
        Lm=laguerre(m)(x)
        print "Ln, Lm",Ln,Lm
        sum=sum+Ln*Lm*glwts[i]
    return sum


def funproj(fun,n,order):
    [glpts,glwts]=la_roots(order,0)
    sum=0.
    for i in range(len(glpts)):
        x=glpts[i]
        Ln=laguerre(n)(x)
        fval=fun(x)*exp(abs(x)/2.)
        sum=sum+fval*Ln*glwts[i]
    return sum

def fundecomp(fun,order):
    glorder=2*order
    retarray=zeros(order)
    for i in range(order):
        retarray[i]=funproj(fun,i,glorder)
    return retarray

def bfrep(fdecomp,x):
    sum=0.
    for i in range(len(fdecomp)):
        sum=sum+fdecomp[i]*laguerre(i)(x)*exp(-abs(x)/2.)
    return sum

def testfun(x):
    #return exp(-pow(x,2))#laguerre(1)(x)*exp(-abs(x)/2.)#exp(-pow(x,2))
    (Ai,Aip,Bi,Bip)=airy(x)
    return Ai

#print twofunctionintegral(2,2)
fdecomp=fundecomp(testfun,10)
xvals=arange(0.,10.,.1)
fvals=testfun(xvals)
f2vals=bfrep(fdecomp,xvals)

plot(xvals,fvals)
plot(xvals,f2vals)
show()
