def compareplots():
    plotwf_sum("tmp/boundstate.dat",-exp(-1j*En*.01),"tmp/wf1.dat",1.)
    plotwf("wf1.dat")
    show()
