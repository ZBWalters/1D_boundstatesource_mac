from numpy import *
from pylab import *
from math import pi
import glob
import matplotlib.cm as cm#color maps
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
#import matplotlib.animation as animation#animation requires animation module
#to be installed on the computer (not true at work, but these functions work
#on my laptop -ZW)


def plotwf(filename,clr='b'):
    array1=genfromtxt(filename)
    (n1,n2)=shape(array1)
    cmplxarray=array1[:,1]+1j*array1[:,2]
    xarray=array1[:,0]
    rangearray=.5*(1+angle(cmplxarray)/pi)
    colorarray=cm.hsv(rangearray,1)
    ##colorplot(xarray,cmplxarray)
    #figure()
    #scatter(xarray,abs(cmplxarray),color=colorarray)
    ##for i in range(n1):
    ##    scatter(xarray[i],abs(cmplxarray[i]),color=colorarray[i])
    #show()

    #clf()
    #print "clr=",clr
    if(clr=='color'):
        #print "color plot"
        colorplot(xarray,cmplxarray)
    else:
        componentplot(xarray,cmplxarray,clr)
    return

def plotwf_times_fun(filename,fn):
    array1=genfromtxt(filename)
    (n1,n2)=shape(array1)
    cmplxarray=array1[:,1]+1j*array1[:,2]
    xarray=array1[:,0]
    fnarray=fn(xarray)
    for i in range(len(cmplxarray)):
        cmplxarray[i]=cmplxarray[i]*fnarray[i]
    rangearray=.5*(1+angle(cmplxarray)/pi)
    colorarray=cm.hsv(rangearray,1)
    ##colorplot(xarray,cmplxarray)
    #figure()
    #scatter(xarray,abs(cmplxarray),color=colorarray)
    ##for i in range(n1):
    ##    scatter(xarray[i],abs(cmplxarray[i]),color=colorarray[i])
    #show()

    #clf()
    

    colorplot(xarray,cmplxarray)
    #componentplot(xarray,cmplxarray,'b')
    return

def plotwf_gauge(statefilename,tplot,F0,En=-.669818,wL=.0565,clr='b'):
    A0=F0/wL
    array1=genfromtxt(statefilename)
    (n1,n2)=shape(array1)
    xarray=array1[:,0]
    cmplxarray=array1[:,1]+1j*array1[:,2]
    for i in range(len(xarray)):
        cmplxarray[i]=cmplxarray[i]*exp(1j*A0*sin(wL*tplot)*xarray[i])*exp(-1j*En*tplot)
    componentplot(xarray,cmplxarray,clr)


def plotwf_minus_gauge(filename,statefilename,tplot,F0,En=-.669818,wL=.0565,clr='b'):
    A0=F0/wL
    array1=genfromtxt(filename)
    (n1,n2)=shape(array1)
    cmplxarray=array1[:,1]+1j*array1[:,2]
    xarray=array1[:,0]

    array2=genfromtxt(statefilename)
    (n1,n2)=shape(array2)
    cmplxarray2=array2[:,1]+1j*array2[:,2]
    
    diffarray=cmplxarray2*0j

    for i in range(len(xarray)):
        cmplxarray2[i]=cmplxarray2[i]*exp(1j*A0*sin(wL*tplot)*xarray[i])*exp(-1j*En*tplot)
        diffarray[i]=cmplxarray[i]-cmplxarray2[i]

    componentplot(xarray,diffarray,clr)

def plotwf_abs_color(filename,lstyle='-',plotcolor='b',lgnd=''):
    array1=genfromtxt(filename)
    (n1,n2)=shape(array1)
    cmplxarray=array1[:,1]+1j*array1[:,2]
    xarray=array1[:,0]
    rangearray=.5*(1+angle(cmplxarray)/pi)
    colorarray=cm.hsv(rangearray,1)
    ##colorplot(xarray,cmplxarray)
    #figure()
    #scatter(xarray,abs(cmplxarray),color=colorarray)
    ##for i in range(n1):
    ##    scatter(xarray[i],abs(cmplxarray[i]),color=colorarray[i])
    #show()

    #clf()
    

    #colorplot(xarray,cmplxarray)
    absplot(xarray,cmplxarray,lstyle,plotcolor,lgnd)
    return

def plotwf_components_color(filename,plotcolor='b',lgnd=''):
    array1=genfromtxt(filename)
    (n1,n2)=shape(array1)
    cmplxarray=array1[:,1]+1j*array1[:,2]
    xarray=array1[:,0]
    rangearray=.5*(1+angle(cmplxarray)/pi)
    colorarray=cm.hsv(rangearray,1)
    ##colorplot(xarray,cmplxarray)
    #figure()
    #scatter(xarray,abs(cmplxarray),color=colorarray)
    ##for i in range(n1):
    ##    scatter(xarray[i],abs(cmplxarray[i]),color=colorarray[i])
    #show()

    #clf()
    

    #colorplot(xarray,cmplxarray)
    componentplot(xarray,cmplxarray,plotcolor,lgnd)
    return


def plotwf_sum(filename1,coeff1,filename2,coeff2,clr='k',lgnd='',):
    array1=genfromtxt(filename1)
    (n1,n2)=shape(array1)
    array2=genfromtxt(filename2)
    
    cmplxarray=coeff1*(array1[:,1]+1j*array1[:,2])+coeff2*(array2[:,1]+1j*array2[:,2])
    xarray=array1[:,0]
    rangearray=.5*(1+angle(cmplxarray)/pi)
    colorarray=cm.hsv(rangearray,1)
    colorplot(xarray,cmplxarray)
    #componentplot(xarray,cmplxarray,clr,lgnd)
    return

def colorplot(x,z,lwidth=2.,alph=1.):
    #tiny=1e-10
    #y=log(abs(z)+tiny)
    y=abs(z)
    theta=angle(0.5*(z[:-1]+z[1:]))
    colorarray=cm.hsv(.5*(1+theta/pi))

    # Create a set of line segments so that we can color them individually
    # This creates the points as a N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line
    # collection needs to be numlines x points per line x 2 (x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    lc = LineCollection(segments, colors=colorarray,linewidth=lwidth,alpha=alph)

    #fig1 = plt.figure()
    plt.gca().add_collection(lc)
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(),y.max())
    #plt.show()

def componentplot(x,z,plotcolor='k',lgnd=''):
    plot(x,real(z),'--',color=plotcolor,linewidth=.5)
    plot(x,imag(z),':',color=plotcolor,linewidth=.5)

def absplot(x,z,lstyle='-',plotcolor='k',lgnd=''):
    plot(x,abs(z),linestyle=lstyle,color=plotcolor,linewidth=.5,label=lgnd)
    
    
    
def init():
    xdat=[0]
    ydat=[0]
    plot(xdat,ydat)

def globplot(globstring,clr='color'):
    figure()
    filelist=glob.glob(globstring)
    natsort(filelist)
    for filename in filelist:
        print "plotting file "+filename
        plotwf(filename,clr)
    show()

def globplot_savefigs(globstring,clr='color'):
    filelist=glob.glob(globstring)
    natsort(filelist)
    fignamestring=""
    for filename in filelist:
        figure()
        print "plotting file "+filename
        plotwf(filename,clr)
        figname=filename[:-4]+".png"
        savefig(figname)
        close()
        fignamestring+=figname+" "
    print "command for animation:"
    print "convert -delay 10 -loop 0 "+fignamestring+"animation.gif"
    

def globplot_lastn(globstring,n=1,clr='color'):
    figure()
    filelist=glob.glob(globstring)
    natsort(filelist)
    for filename in filelist[-n:]:
        print "showing file "+filename
        plotwf(filename,clr)
    show()

def plotlist(filelist):
    figure()
    for i in range(len(filelist)):
        filename=filelist[i]
        print "showing file "+filename
        plotwf(filename)
    show()


# ---------------------------------------------------------
# natsort.py: Natural string sorting.
# ---------------------------------------------------------

# By Seo Sanghyeon.  Some changes by Connelly Barnes.

def try_int(s):
    "Convert to integer if possible."
    try: return int(s)
    except: return s

def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
    "Natural string comparison, ignores case."
    return natcmp(a.lower(), b.lower())

def natsort(seq, cmp=natcmp):
    "In-place natural string sort."
    seq.sort(cmp)
    
def natsorted(seq, cmp=natcmp):
    "Returns a copy of seq, sorted by natural string sort."
    import copy
    temp = copy.copy(seq)
    natsort(temp, cmp)
    return temp


###############functions for animation

def animate(i,filelist):
    filename=filelist[i];
    plotwf(filename)

def globanimate(globstring):
    fig=plt.figure()
    window=fig.add_subplot(1,1,1)
    filelist=glob.glob(globstring)
    natsort(filelist)
    ani=animation.FuncAnimation(fig,animate,len(filelist),fargs=[filelist])
    show()

################################functions for a particular comparison
def compareplots():
    En=-.669
    plotwf_sum("tmp/boundstate.dat",-exp(-1j*En*.001),"tmp/wf1.dat",1.)
    plotwf("wf1.dat")
    show()