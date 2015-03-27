from numpy import *
from pylab import *

def adjustphase(mat):
    nmat=len(mat)
    retmat=zeros(len(mat))*0j;
    for i in range(len(mat)):
        retmat[i]=mat[i]*((-1.)**i)
    return retmat

def recenter(mat):
    nmat=len(mat)
    retmat=zeros(len(mat))*0j;
    for i in range(len(mat)):
        retmat[(i+nmat/2)%nmat]=mat[i]
    return retmat

def postprocess(mat):
    nmat=len(mat)
    retmat=zeros(len(mat))*0j;
    for i in range(len(mat)):
        retmat[(i+nmat/2)%nmat]=mat[i]*((-1.)**i)/nmat
    return retmat

def kpts(nxpts):
    return arange(nxpts)-nxpts/2

def xpts(xmin,xmax,nxpts):
    dx=(xmax-xmin)/(nxpts-1.)
    return arange(xmin,xmax+dx,dx)

def func(x):
    return x*exp(-abs(x))*(1.+0j)

def funcvals(xmin,xmax,nxpts):
    xarray=xpts(xmin,xmax,nxpts)
    return func(xarray)

array128=genfromtxt("wf_grid.128")
array256=genfromtxt("wf_grid.256")

k128=kpts(128)
x128=xpts(-10,10.,128)#array128[:,0]
psi128=funcvals(-10.,10.,128)#array128[:,1]+1j*array128[:,2]

k256=kpts(256)
x256=xpts(-10.,10.,256)#array256[:,0]
psi256=funcvals(-10.,10.,256)#array256[:,1]+1j*array256[:,2]

psi102=funcvals(-10.,10.,102)
x102=xpts(-10.,10.,101)
k102=kpts(102)

fft102=fft(psi102)
fft128=fft(psi128)
fft256=fft(psi256)


#plot(x102,real(postprocess(fft102)))
#plot(x128,real(postprocess(fft128)))
#plot(x256,real(postprocess(fft256)))
#show()

plot(k102,imag(postprocess(fft102)))
plot(k128,imag(postprocess(fft128)))
plot(k256,imag(postprocess(fft256)))
show()
