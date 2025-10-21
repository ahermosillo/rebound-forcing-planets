#!/usr/bin/python
import sys, getopt
import numpy as np
#import matplotlib.pyplot as plt
from os.path import exists



def mod2pi(x):
    '''
    input:
        x = any angle in radians
    output
        an angle in radians re-centered from 0-2pi
    '''
    nmax=len(x)
    for n in range(0,nmax):
        while(x[n]>2.*np.pi):
            x[n]+=-2.*np.pi
        while(x[n]<0.):
            x[n]+=2.*np.pi
    return x

prefix = 'lowe-164-'
pl = 4 #Neptune's particle id number
tmin = 0. #1.99e8

aei = open("all-end-states.txt","a")
aei.write("\#20 Myr average a(1), e(2), i(3), final osculating a(4), e(5), i(6)\n")
taei = open("2-1-end-states.txt","a")
taei.write("\#20 Myr average a(1), e(2), i(3), final osculating a(4), e(5), i(6), delta phi(7), phibar(8), number of librating windows(9), avg. window dphi(10)\n")


for diri in range(101,102):

    infile = prefix + str(diri) + '/high-res-end-orbits.txt'
    print(infile)
    if( not (exists(infile))):
        continue
    data = np.genfromtxt(infile, names=['id', 't', 'a', 'e', 'inc', 'node', 'peri', 'MA'])

    time = data['t'][np.logical_and(data['id']==pl,data['t'] > tmin )]
    
    lpl = (data['node'][np.logical_and(data['id']==pl,data['t'] > tmin )] + 
        data['peri'][np.logical_and(data['id']==pl,data['t'] > tmin )] + 
        data['MA'][np.logical_and(data['id']==pl,data['t'] > tmin )])
    lpl = mod2pi(lpl)
    #wtpl = data['node'][(data['id']==pl)] + data['peri'][(data['id']==pl)]
    #nodepl =  data['node'][(data['id']==pl)]

    anep = np.mean(data['a'][np.logical_and(data['id']==pl,data['t'] > tmin )])
    ares = anep*(np.power(2.,(2./3.)))
    print(ares)
    

    for tp in range(101,251):
        atp = data['a'][np.logical_and(data['id']==tp, data['t'] > tmin )]
        etp = data['e'][np.logical_and(data['id']==tp, data['t'] > tmin )]
        itp = data['inc'][np.logical_and(data['id']==tp, data['t'] > tmin )]
        ltp = (data['node'][np.logical_and(data['id']==tp, data['t'] > tmin )] + 
                data['peri'][np.logical_and(data['id']==tp, data['t'] > tmin )] + 
                data['MA'][np.logical_and(data['id']==tp, data['t'] > tmin )])
        ltp = mod2pi(ltp)
        wttp = (data['node'][np.logical_and(data['id']==tp, data['t'] > tmin )] + 
                data['peri'][np.logical_and(data['id']==tp, data['t'] > tmin )])
        wttp = mod2pi(wttp)
        #nodetp = data['node'][np.logical_and(data['id']==tp, data['t'] > tmin )]


        phi = 2.*ltp - lpl - wttp
        phi = mod2pi(phi)

        nmax = len(atp)
        print("tp %d nmax %d" %(tp,nmax))
        if(nmax >100):

            abar = np.mean(atp)
            ebar = np.mean(etp)
            ibar = np.mean(itp)*180./np.pi

            aei.write("%d %10f %10f %10f %10f %10f %10f\n" % (tp, abar, ebar, ibar, atp[-1], etp[-1], itp[-1]))

            da = np.abs(abar-ares)

            if(da < 0.3):
                nmax=len(phi)

                phi = phi*180/np.pi
                nwin = int(np.floor(nmax/20.))
                lib = 0
                dphibar = 0.
                for n in range(0,nmax+1-nwin,nwin):
                    n2 = n+nwin
                    phit = phi[n:n2]
                    dphit = 0.5*(np.amax(phit) - np.amin(phit))
                    if(dphit < 172):
                        lib+=1
                        dphibar+=dphit
    
                phimin = np.amin(phi)
                phimax = np.amax(phi)
    
                dphi = (phimax - phimin)/2.
                phibar = np.mean(phi)
    
                if(lib>10):
                    dphibar = dphibar/lib
                    taei.write("%d %10f %10f %10f %10f %10f %10f %10f %10f %d %f\n" % 
                            (tp, abar, ebar, ibar, atp[-1], etp[-1], itp[-1], dphi, phibar, lib, dphibar))

      
aei.close()
taei.close()

quit()
