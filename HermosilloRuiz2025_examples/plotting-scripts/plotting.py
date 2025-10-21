import matplotlib.pyplot as plt
import numpy as np
import rebound
# import pyorb
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# from numba import jit
import glob
plt.style.use('paper_labels_colors.mplstyle')

def get_data(filename,nP, rebound = True,names = ["iP", "Time", "a", "e", "i","peri", "node", "f"],tmin=None,tmax=None):
    if rebound:
        data = np.genfromtxt(filename,names=names)
        for i in range(1,nP+1):
            pl = data[data['iP'] == i]
            if (tmin==None) & (tmax==None):
                pass
            else:
                pl = pl[(pl['Time']>=tmin) & (pl['Time']<=tmax)]
            if i==1:
                planets = [pl]
            else:
                planets = np.append(planets,[pl],axis=0)
        return planets
    else:
        for i in range(nP):
            aeifilename = "{}".format(filename[i])           
            # get rid of extra header field that messes things up
            with open(aeifilename,'r',encoding='utf-8') as file:
                tempdata = file.readlines()
            if(tempdata[3] == "    Time (years)        a        e        i      peri     node       f        M        mass     \n"):
                tempdata[3] = "    Time                a        e        i      peri     node       f        M        mass     \n"
            with open(aeifilename, 'w', encoding='utf-8') as file:
                file.writelines(tempdata)                                
            pl = np.genfromtxt(aeifilename,skip_header=3,names=True)
#             print(pl)
            if (tmin==None) & (tmax==None):
                pass
            else:
                pl = pl[(pl['Time']>=tmin) & (pl['Time']<=tmax)]
            if i==0:
                planets = [pl]
            else:
                planets = np.append(planets,[pl],axis=0)
        return np.array(planets)


def exponential(data,param,delta,tau,trange):
    t_dat = data['Time']
    if trange is None:
        t = t_dat
    else:
        t = t_dat[(t_dat>trange[0])&(t_dat<trange[1])]
    eq = (data[param][0] + delta) - delta*np.exp(-t/tau)
    return t, eq

def linear(data,param,delta,tau,trange):
    t_dat = data['Time']
    if trange is None:
        t = t_dat
        print("trange was None")
    else:
        print("trange was not None")
        t = t_dat[(t_dat>trange[0])&(t_dat<trange[1])]
    eq = (data[param][0]) + delta*(t/tau)
    return t, eq

def sinusoidal(data,param,delta,tau,trange):
    t_dat = data['Time']
    if trange is None:
        t = t_dat
    else:
        t = t_dat[(t_dat>trange[0])&(t_dat<trange[1])]
#     print(t)
    eq = data[param][0] + delta*np.sin(2.0*np.pi*t/tau)
    return t, eq

def logarithmic(data,param,delta,tau,trange):
    t_dat = data['Time']
    if trange is None:
        t = t_dat
    else:
        t = t_dat[(t_dat>trange[0])&(t_dat<trange[1])]
    eq = (data[param][0]) + delta*np.log(t/tau + 1)
    return t, eq

def get_equation(data,param,delta,tau,eqBool,trange):
    t_dat = data['Time']
    if trange is None:
        t = t_dat
    else:
        t = t_dat[(t_dat>trange[0])&(t_dat<trange[1])]
    if eqBool==0:
        return exponential(data,param,delta,tau,trange)
    elif eqBool==1:
        return linear(data,param,delta,tau,trange)
    elif eqBool==2:
#         print("sinusoidal function")
        return sinusoidal(data,param,delta,tau,trange)
    elif eqBool==3:
        return logarithmic(data,param,delta,tau,trange)
    else:
#         print("not a valid equation")
        # return empty lists to plot no curve
        # return [],[]
        return t, np.ones(len(t))*np.nan

def calc_xyz_plane(a, e, i, Om, w, f):
    """calculating x y and z given 6 orbital parameters"""
    P1 = np.array([np.cos(w), -np.sin(w), 0, np.sin(w), np.cos(w), 0, 0, 0, 1]).reshape(3,3)
    P2 = np.array([1, 0, 0, 0, np.cos(i), -np.sin(i), 0, np.sin(i), np.cos(i)]).reshape(3,3)
    P3 = np.array([np.cos(Om), -np.sin(Om), 0, np.sin(Om), np.cos(Om), 0, 0, 0, 1]).reshape(3,3)
    r = a*(1-e*e)/(1 + e*np.cos(f))
    rf = np.zeros(3).reshape(3,1)
    rf[0] = r*np.cos(f)
    rf[1] = r*np.sin(f)
    
    t1 = np.matmul(P3,P2)
    t2 = np.matmul(t1, P1)
    t3 = np.matmul(t2, rf)
    
    return t3


def publication_plots_ae(data,indxPlanets,finalParam,tau,expEq,force_trange=None,rebound=True,verbose=True):
    print(f"plot only {indxPlanets} planets")
    nP = len(indxPlanets)
    fig,axs = plt.subplots(nP,2,figsize = (4,3),sharex=True)
    axs = axs.ravel()

    das = np.zeros(nP)
    des = np.zeros(nP)
    # print(das)
    # print(len(das))

    t_dat = data[0]['Time']
    if force_trange is not None: 
        t_dat = t_dat[(t_dat>force_trange[0])&(t_dat<force_trange[1])]
    else:
        pass
    
    nt = len(t_dat)
    aeqs = np.empty((nP,nt))
    eeqs = np.empty((nP,nt))
    
    for j,i in enumerate(indxPlanets):
        # print("j",j)
        # print("i",i)
        das[j] = finalParam[0][i] - data[i]['a'][0]
        des[j] = finalParam[1][i] - data[i]['e'][0]
        # print(data[i])
        # print(das[j])
        # print(tau[i])
        # print(expEq[i])
        aeq = get_equation(data[i],'a',das[j],tau[i][0],expEq[i][0],force_trange)
        eeq = get_equation(data[i],'e',des[j],tau[i][1],expEq[i][1],force_trange)
        
        t_eq,aeqs[j] =aeq[0],np.array(aeq[1])
        t_eq,eeqs[j] = eeq[0],np.array(eeq[1])

        if verbose:
            print("for planet {}: a0 = {}; e0 = {}; i0 = {};".format(i,data[i]['a'][0],
                    data[i]['e'][0],data[i]['i'][0]))
            print(f"w0 = {data[i]['peri'][0]*180/np.pi}; Om0 = {data[i]['node'][0]*180/np.pi}; f0 = {data[i]['f'][0]*180/np.pi}")
            print("da = {}; de = {};".format(das[j],des[j]))
            print(f"taua = {tau[i][0]/(2*np.pi)} yr; taue = {tau[i][1]/(2*np.pi)} yr, taui = {tau[i][2]/(2*np.pi)} yr, tauw = {tau[i][3]/(2*np.pi)} yr; tauOm = {tau[i][4]/(2*np.pi)} yr")
            
    if rebound: 
        t_dat /= (2*np.pi)
        t_eq/=(2*np.pi)
    t_dat/=1e6
    t_eq/=1e6
    
    for i,j in enumerate(indxPlanets):
        simulationcolor = '#EE6C4D'
        dampfunctioncolor = '#750706'
        linestyle = '-'

        if j==2:
            dampfunctioncolor = '#354E6E'
            linestyle = '--'

        axs[i*2].plot(t_dat,data[j]['a'],marker='.',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        axs[i*2].plot(t_eq,aeqs[i],color=dampfunctioncolor,linestyle=linestyle,rasterized=True)

        axs[i*2+1].plot(t_dat,data[j]['e'],marker='.',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        axs[i*2+1].plot(t_eq,eeqs[i],color=dampfunctioncolor,linestyle=linestyle,rasterized=True)
        
        # if not np.isnan(max(data[j]['e'])):
        #     axs[i*3+1].set_ylim(0,max(data[j]['e'])+0.1*max(data[j]['e']))
        # axs[i*3+1].set_ylim(0,0.11)

        # if not np.isnan(max(data[j]['i']*180/np.pi)):
        #     axs[i*3+2].set_ylim(0,max(data[j]['i']*180/np.pi)+0.1*max(data[j]['i']*180/np.pi))
        # axs[i*3+2].set_ylim(-1,12)
        
        # axs[i*2].yaxis.set_minor_locator(AutoMinorLocator(2))
        # axs[i*2+1].yaxis.set_minor_locator(AutoMinorLocator(5))
        
        axs[i*2].set_ylabel("$a$ (au)")
        axs[i*2+1].set_ylabel("$e$")
  
    return fig, axs

def publication_plots(data,indxPlanets,finalParam,tau,expEq,force_trange=None,rebound=True,verbose=True):
    print(f"plot only {indxPlanets} planets")
    nP = len(indxPlanets)
    fig,axs = plt.subplots(nP,3,figsize = (6,1.5*nP),sharex=True)
    axs = axs.ravel()

    das = np.zeros(nP)
    des = np.zeros(nP)
    dis = np.zeros(nP)

    t_dat = data[0]['Time']
    if force_trange is not None: 
        t_dat = t_dat[(t_dat>force_trange[0])&(t_dat<force_trange[1])]
    else:
        pass
    
    nt = len(t_dat)
    aeqs = np.empty((nP,nt))
    eeqs = np.empty((nP,nt))
    ieqs = np.empty((nP,nt))
    
    for j,i in enumerate(indxPlanets):
    #         print(i)
        das[j] = finalParam[0][i] - data[i]['a'][0]
        des[j] = finalParam[1][i] - data[i]['e'][0]
        dis[j] = finalParam[2][i] - data[i]['i'][0]

        aeq = get_equation(data[i],'a',das[j],tau[i][0],expEq[i][0],force_trange)
        eeq = get_equation(data[i],'e',des[j],tau[i][1],expEq[i][1],force_trange)
        ieq = get_equation(data[i],'i',dis[j],tau[i][2],expEq[i][2],force_trange)
        t_eq,aeqs[j] =aeq[0],np.array(aeq[1])
        t_eq,eeqs[j] = eeq[0],np.array(eeq[1])
        t_eq,ieqs[j] = ieq[0],np.array(ieq[1])%(2*np.pi)

        if verbose:
            print("for planet {}: a0 = {}; e0 = {}; i0 = {};".format(i,data[i]['a'][0],
                    data[i]['e'][0],data[i]['i'][0]))
            print(f"w0 = {data[i]['peri'][0]*180/np.pi}; Om0 = {data[i]['node'][0]*180/np.pi}; f0 = {data[i]['f'][0]*180/np.pi}")
            print("da = {}; de = {}; di = {};".format(das[j],des[j],dis[j]))
            print(f"taua = {tau[i][0]/(2*np.pi)} yr; taue = {tau[i][1]/(2*np.pi)} yr, taui = {tau[i][2]/(2*np.pi)} yr, tauw = {tau[i][3]/(2*np.pi)} yr; tauOm = {tau[i][4]/(2*np.pi)} yr")
            
    if rebound: 
        t_dat /= (2*np.pi)
        t_eq/=(2*np.pi)
    t_dat/=1e6
    t_eq/=1e6

    for i,j in enumerate(indxPlanets):
        simulationcolor = '#EE6C4D'
        dampfunctioncolor = '#750706'
        incdat = data[j]['i']

        if j==2:
            dampfunctioncolor = '#354E6E'

        axs[j*3].plot(t_dat,data[j]['a'],marker='.',color=simulationcolor,markersize=0.5,rasterized=True)
        axs[j*3].plot(t_eq,aeqs[j],color=dampfunctioncolor,linestyle='--',markersize=3,rasterized=True)
        axs[j*3+1].plot(t_dat,data[j]['e'],marker='.',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        axs[j*3+1].plot(t_eq,eeqs[j],color=dampfunctioncolor,linestyle='--',markersize=3,rasterized=True)
        # if not np.isnan(max(data[j]['e'])):
        #     axs[i*3+1].set_ylim(0,max(data[j]['e'])+0.1*max(data[j]['e']))
        # axs[i*3+1].set_ylim(0,0.11)

        axs[j*3+2].plot(t_dat,incdat*180/np.pi,marker='.',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        axs[j*3+2].plot(t_eq,ieqs[j]*180/np.pi,color=dampfunctioncolor,linestyle='--',markersize=3,rasterized=True)
        # if not np.isnan(max(data[j]['i']*180/np.pi)):
        #     axs[i*3+2].set_ylim(0,max(data[j]['i']*180/np.pi)+0.1*max(data[j]['i']*180/np.pi))
        # axs[i*3+2].set_ylim(-1,12)
        
        axs[i*3].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*3+1].yaxis.set_minor_locator(AutoMinorLocator(5))
        axs[i*3+2].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*3].xaxis.set_minor_locator(AutoMinorLocator(3))
        axs[i*3].xaxis.set_major_locator(MultipleLocator(30))

        axs[i*3].set_ylabel("$a$ (au)")
        axs[i*3+1].set_ylabel("$e$")
        axs[i*3+2].set_ylabel("$i$ (deg)")
    return fig, axs

def plot_planets(data,finalParam,tau,expEq,force_trange=None,rebound=True,verbose=True):
    nP = len(data)
    print(f"there are {nP} planets")
    fig,axs = plt.subplots(nP*2,3,figsize=(6,nP*3),sharex=True,layout='constrained')
    axs = axs.ravel()
    
    das = np.zeros(nP)
    des = np.zeros(nP)
    dis = np.zeros(nP)
    dOms = np.zeros(nP)
    dws = np.zeros(nP)

    t_dat = data[0]['Time']
    if force_trange is not None: 
        t = t_dat[(t_dat>force_trange[0])&(t_dat<force_trange[1])]
    else:
        t = t_dat

    nt = len(t)
    aeqs = np.empty((nP,nt))
    eeqs = np.empty((nP,nt))
    ieqs = np.empty((nP,nt))
    Omeqs = np.empty((nP,nt))
    weqs = np.empty((nP,nt))
    print(aeqs.shape)

    simulationcolor = '#EE6C4D'
    dampfunctioncolor = '#750706'
    
    for i in range(nP):
#         print(i)
        das[i] = finalParam[0][i] - data[i]['a'][0]
        des[i] = finalParam[1][i] - data[i]['e'][0]
        dis[i] = finalParam[2][i] - data[i]['i'][0]
        dws[i] = finalParam[3][i] - data[i]['peri'][0]
        dOms[i] = finalParam[4][i] - data[i]['node'][0]

        aeq = get_equation(data[i],'a',das[i],tau[i][0],expEq[i][0],force_trange)
        eeq = get_equation(data[i],'e',des[i],tau[i][1],expEq[i][1],force_trange)
        ieq = get_equation(data[i],'i',dis[i],tau[i][2],expEq[i][2],force_trange)
        weq = get_equation(data[i],'peri',dws[i],tau[i][3],expEq[i][3],force_trange)
        Omeq = get_equation(data[i],'node',dOms[i],tau[i][4],expEq[i][4],force_trange)

        t_eq,aeqs[i] =aeq[0],np.array(aeq[1])
        t_eq,eeqs[i] = eeq[0],np.array(eeq[1])
        t_eq,ieqs[i] = ieq[0],np.array(ieq[1])%(2*np.pi)
        t_eq,weqs[i] = weq[0],np.array(weq[1])%(2*np.pi)
        t_eq,Omeqs[i] = Omeq[0],np.array(Omeq[1])%(2*np.pi)

        if verbose:
            print("for planet {}: a0 = {}; e0 = {}; i0 = {};".format(i,data[i]['a'][0],
                                            data[i]['e'][0],data[i]['i'][0]))
            print(f"w0 = {data[i]['peri'][0]*180/np.pi}; Om0 = {data[i]['node'][0]*180/np.pi}; f0 = {data[i]['f'][0]*180/np.pi}")
            print("da = {}; de = {}; di = {};".format(das[i],des[i],dis[i]))
            print(f"taua = {tau[i][0]/(2*np.pi)} yr; taue = {tau[i][1]/(2*np.pi)} yr, taui = {tau[i][2]/(2*np.pi)} yr, tauw = {tau[i][3]/(2*np.pi)} yr; tauOm = {tau[i][4]/(2*np.pi)} yr")
        
    if rebound: 
        t_dat /= (2*np.pi)
        t_eq/=(2*np.pi)
    t_dat/=1e6
    t_eq/=1e6
    for j in range(nP):
        print(j)
        incdat = data[j]['i']%(2*np.pi)
        peridat = data[j]['peri']%(2*np.pi)
        nodedat = data[j]['node']%(2*np.pi)
        Omeq = Omeqs[j]%(2*np.pi)
 
#         if any(abs(np.diff(incdat))>350):
#             incdat[incdat>np.pi] = (incdat[incdat>np.pi]+np.pi)%(-np.pi)
#         else:
#             pass
        # if any(abs(np.diff(peridat))>350*np.pi/180):
        #     peridat[peridat>np.pi] = (peridat[peridat>np.pi]+np.pi)%(-np.pi)
        # else:
        #     pass
        # print(abs(np.diff(nodedat)))
        if any(abs(np.diff(nodedat))>350*np.pi/180):
            # print("true")
            nodedat[nodedat>np.pi] = (nodedat[nodedat>np.pi]-2*np.pi)
        else:
            pass
        if any(abs(np.diff(Omeq))>350*np.pi/180):
            print("true")
            Omeq[Omeq>np.pi] = (Omeq[Omeq>np.pi]-2*np.pi)
        else:
            pass
#         print(j)
        if force_trange is not None:
            t_i = force_trange[0]
            t_f = force_trange[1]
        else:
            pass
        print(t_dat)
        print(data[j]['a'])
        axs[j*6].plot(t_dat,data[j]['a'],marker='.',color=simulationcolor,markersize=0.5,rasterized=True)
        axs[j*6].plot(t_eq,aeqs[j],color=dampfunctioncolor,linestyle='--',markersize=3,rasterized=True)
        axs[j*6+1].plot(t_dat,data[j]['e'],marker='.',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        axs[j*6+1].plot(t_eq,eeqs[j],color=dampfunctioncolor,linestyle='--',markersize=3,rasterized=True)
        if not np.isnan(max(data[j]['e'])):
            axs[j*6+1].set_ylim(min(data[j]['e'])-0.1*min(data[j]['e']),max(data[j]['e'])+0.1*max(data[j]['e']))
        axs[j*6+2].plot(t_dat,incdat*180/np.pi,marker='.',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        axs[j*6+2].plot(t_eq,ieqs[j]*180/np.pi,color=dampfunctioncolor,linestyle='--',markersize=3,rasterized=True)
        if not np.isnan(max(data[j]['i']*180/np.pi)):
            axs[j*6+2].set_ylim(min(data[j]['i']*180/np.pi)-0.1*min(data[j]['i']*180/np.pi),max(data[j]['i']*180/np.pi)+0.1*max(data[j]['i']*180/np.pi))
        axs[j*6+3].plot(t_dat,peridat*180/np.pi,marker='o',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        axs[j*6+4].plot(t_dat,nodedat*180/np.pi,marker='o',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        axs[j*6+5].plot(t_dat,data[j]['f']*180/np.pi,marker='o',color=simulationcolor,linestyle='',markersize=0.5,rasterized=True)
        
        axs[j*6+3].plot(t_eq,weqs[j]*180/np.pi,color=dampfunctioncolor,marker="o",linestyle='',markersize=0.1,rasterized=True)
        axs[j*6+4].plot(t_eq,Omeq*180/np.pi,color=dampfunctioncolor,marker="o",linestyle='',markersize=0.1,rasterized=True)
        #axs[j*3+5].plot(t,[j],marker='.',color=dampfunctioncolor,linestyle='',markersize=3)
        axs[j*6+1].set_ylabel("$e$")
        axs[j*6+2].set_ylabel("$i$ (deg)")
        axs[j*6].set_ylabel("$a$ (au)")
        axs[j*6+3].set_ylabel("$\omega$ (deg)")
        axs[j*6+4].set_ylabel("$\Omega$ (deg)")
        axs[j*6+5].set_ylabel("$f$ (deg)")
        axs[i*6].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+1].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+2].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+3].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+4].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+5].yaxis.set_minor_locator(AutoMinorLocator(2))
        # axs[i*6].xaxis.set_minor_locator(AutoMinorLocator(3))
        # axs[i*6].xaxis.set_major_locator(MultipleLocator(15))
    return fig, axs

def plot_vusr_novusr_ex(file1,file2,finalParam,tau,expEq,force_trange=None,rebound=True):
    no_vusr_data = file1
    data=file2
    nt = len(data[0])
    nP = len(data)
    print(f"there are {nP} planets")
    fig,axs = plt.subplots(nP*2,3,figsize=(6,nP*3),sharex=True,layout='constrained')
    axs = axs.ravel()
    
    das = np.zeros(nP)
    des = np.zeros(nP)
    dis = np.zeros(nP)
    dOms = np.zeros(nP)
    dws = np.zeros(nP)

    t_dat = data[0]['Time']
    if force_trange is not None: 
        t = t_dat[(t_dat>force_trange[0])&(t_dat<force_trange[1])]
    else:
        t = t_dat

    nt = len(t)
    aeqs = np.empty((nP,nt))
    eeqs = np.empty((nP,nt))
    ieqs = np.empty((nP,nt))
    Omeqs = np.empty((nP,nt))
    weqs = np.empty((nP,nt))
    # print(aeqs.shape)

    simulationcolor = '#EE6C4D'
    dampfunctioncolor = '#750706'
    
    for i in range(nP):
#         print(i)
        das[i] = finalParam[0][i] - data[i]['a'][0]
        des[i] = finalParam[1][i] - data[i]['e'][0]
        dis[i] = finalParam[2][i] - data[i]['i'][0]
        dws[i] = finalParam[3][i] - data[i]['peri'][0]
        dOms[i] = finalParam[4][i] - data[i]['node'][0]

        aeq = get_equation(data[i],'a',das[i],tau[i][0],expEq[i][0],force_trange)
        eeq = get_equation(data[i],'e',des[i],tau[i][1],expEq[i][1],force_trange)
        ieq = get_equation(data[i],'i',dis[i],tau[i][2],expEq[i][2],force_trange)
        weq = get_equation(data[i],'peri',dws[i],tau[i][3],expEq[i][3],force_trange)
        Omeq = get_equation(data[i],'node',dOms[i],tau[i][4],expEq[i][4],force_trange)

        t_eq,aeqs[i] =aeq[0],np.array(aeq[1])
        t_eq,eeqs[i] = eeq[0],np.array(eeq[1])
        t_eq,ieqs[i] = ieq[0],np.array(ieq[1])%(2*np.pi)
        t_eq,weqs[i] = weq[0],np.array(weq[1])%(2*np.pi)
        t_eq,Omeqs[i] = Omeq[0],np.array(Omeq[1])%(2*np.pi)
    # print(t_dat)
    # print(t_eq)
    # print(t_dat/(2*np.pi))

    if rebound: 
        t_dat /= ((2*np.pi))
        # t_eq/=((2*np.pi))
    # print(t_dat)
    # print(t_eq)
    t_dat/=1e6
    # t_eq/=1e6
    # print(t_dat)
    # print(t_eq)
    for j in range(nP):
        # print(j)
        incdat = data[j]['i']%(2*np.pi)
        peridat = data[j]['peri']%(2*np.pi)
        nodedat = data[j]['node']%(2*np.pi)
        Omeq = Omeqs[j]%(2*np.pi)
 
#         if any(abs(np.diff(incdat))>350):
#             incdat[incdat>np.pi] = (incdat[incdat>np.pi]+np.pi)%(-np.pi)
#         else:
#             pass
        # if any(abs(np.diff(peridat))>350*np.pi/180):
        #     peridat[peridat>np.pi] = (peridat[peridat>np.pi]+np.pi)%(-np.pi)
        # else:
        #     pass
        # print(abs(np.diff(nodedat)))
        if any(abs(np.diff(nodedat))>350*np.pi/180):
            # print("true")
            nodedat[nodedat>np.pi] = (nodedat[nodedat>np.pi]-2*np.pi)
        else:
            pass
        if any(abs(np.diff(Omeq))>350*np.pi/180):
            # print("true")
            Omeq[Omeq>np.pi] = (Omeq[Omeq>np.pi]-2*np.pi)
        else:
            pass
#         print(j)
        if force_trange is not None:
            t_i = force_trange[0]
            t_f = force_trange[1]
        else:
            pass
        # print(t_dat)
        # print(data[j]['a'])
        axs[j*6].scatter(t_dat,data[j]['a'],marker='.',color=simulationcolor,s=0.06)
        axs[j*6].plot(t_eq,aeqs[j],color=dampfunctioncolor,linestyle='--',markersize=3)
        axs[j*6].scatter(t_dat,no_vusr_data[j]['a'],marker=',',s=1,color = '#335471')
        
        axs[j*6+1].scatter(t_dat,data[j]['e'],marker=',',color=simulationcolor,s=0.06)
        axs[j*6+1].scatter(t_dat,no_vusr_data[j]['e'],marker=',',s=2,color = '#335471')
        axs[j*6+1].plot(t_eq,eeqs[j],color=dampfunctioncolor,linestyle='--',markersize=3)
        if not np.isnan(max(data[j]['e'])):
            axs[j*6+1].set_ylim(min(data[j]['e'])-0.1*min(data[j]['e']),max(data[j]['e'])+0.1*max(data[j]['e']))
        
        
        axs[j*6+2].scatter(t_dat,no_vusr_data[j]['i']%(2*np.pi)*180/np.pi,marker=',',s=1,color = '#335471')
        axs[j*6+2].scatter(t_dat,incdat*180/np.pi,marker=',',color=simulationcolor,s=0.06)
        axs[j*6+2].plot(t_eq,ieqs[j]*180/np.pi,color=dampfunctioncolor,linestyle='--',markersize=3)
        if not np.isnan(max(data[j]['i']*180/np.pi)):
            axs[j*6+2].set_ylim(min(data[j]['i']*180/np.pi)-0.1*min(data[j]['i']*180/np.pi),max(data[j]['i']*180/np.pi)+0.1*max(data[j]['i']*180/np.pi))
        axs[j*6+3].scatter(t_dat,peridat*180/np.pi,marker=',',color=simulationcolor,s=0.06)
        axs[j*6+4].scatter(t_dat,no_vusr_data[j]['node']%(2*np.pi)*180/np.pi,marker=',',s=1,color = '#335471')
        axs[j*6+4].scatter(t_dat,nodedat*180/np.pi,marker=',',color=simulationcolor,s=0.06)
        # axs[j*6+5].plot(t_dat,data[j]['f']*180/np.pi,linestyle='-',color=simulationcolor,linewidth=0.1,alpha=0.5)
        axs[j*6+5].scatter(t_dat,data[j]['f']*180/np.pi,marker=',',color=simulationcolor,s=0.06)
        axs[j*6+3].scatter(t_dat,no_vusr_data[j]['peri']%(2*np.pi)*180/np.pi,marker=',',s=1,color = '#335471')
        
        # axs[j*6+5].plot(t_dat,no_vusr_data[j]['f']%(2*np.pi)*180/np.pi,linestyle='-',alpha=0.5,linewidth=0.1)
        axs[j*6+5].scatter(t_dat,no_vusr_data[j]['f']%(2*np.pi)*180/np.pi,marker=',',s=0.06,alpha=0.6,color = '#335471')
        axs[j*6+3].plot(t_eq,weqs[j]*180/np.pi,color=dampfunctioncolor,marker=",",linestyle='',markersize=0.1)
        axs[j*6+4].plot(t_eq,Omeq*180/np.pi,color=dampfunctioncolor,marker=",",linestyle='',markersize=0.1)
        #axs[j*3+5].plot(t,[j],marker='.',color=dampfunctioncolor,linestyle='',markersize=3)
        axs[j*6+1].set_ylabel("$e$")
        axs[j*6+2].set_ylabel("$i$ (deg)")
        axs[j*6].set_ylabel("$a$ (au)")
        axs[j*6+3].set_ylabel("$\omega$ (deg)")
        axs[j*6+4].set_ylabel("$\Omega$ (deg)")
        axs[j*6+5].set_ylabel("$f$ (deg)")
        axs[i*6].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+1].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+2].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+3].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+4].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6+5].yaxis.set_minor_locator(AutoMinorLocator(2))
        axs[i*6].xaxis.set_minor_locator(AutoMinorLocator(3))
        axs[i*6].xaxis.set_major_locator(MultipleLocator(15))

    return fig, axs

def plot_w_damp(file1,file2,indxPlanets,finalParam,tau,expEq,force_trange=None):
    nP = len(indxPlanets)
    fig,axs = plt.subplots(1,figsize = (3.3,3))
    noevData = file1
    data=file2
    nt = len(data[0])

    dOms = np.zeros(nP)
    dws = np.zeros(nP)
    Omeqs = np.empty((nP,nt))
    weqs = np.empty((nP,nt))

    t_dat = data[0]['Time']
    if force_trange is not None: 
        t = t_dat[(t_dat>force_trange[0])&(t_dat<force_trange[1])]
    else:
        t = t_dat

    simulationcolor = '#EE6C4D'
    dampfunctioncolor = '#750706'
    
    for i in range(nP):
#         print(i)
        dws[i] = finalParam[3][i] - data[i]['peri'][0]
        dOms[i] = finalParam[4][i] - data[i]['node'][0]
        weq = np.array(get_equation(data[i],'peri',dws[i],tau[i][3],eqBool=expEq[i][3],trange=force_trange))%(2*np.pi)
        Omeq = np.array(get_equation(data[i],'node',dOms[i],tau[i][4],eqBool=expEq[i][4],trange=force_trange))%(2*np.pi)

        t_eq,weqs[i] = weq[0],np.array(weq[1])%(2*np.pi)
        t_eq,Omeqs[i] = Omeq[0],np.array(Omeq[1])%(2*np.pi)
        
    for j in range(nP):
        # print(j)
        if rebound: 
            t = data[j]['Time']/(2*np.pi)
        else: 
            t = data[j]['Time']
        incdat = data[j]['i']%(2*np.pi)
        peridat = data[j]['peri']%(2*np.pi)
        nodedat = data[j]['node']%(2*np.pi)
        Omeq = Omeqs[j]%(2*np.pi)

    tnoev = noevData[0]['Time']/(2*np.pi)

    axs.plot(t/1e6,weqs[0]*180/np.pi,'.',color=dampfunctioncolor,rasterized=True)
    axs.plot(t/1e6,peridat*180/np.pi,'.',alpha=0.5,rasterized=True,color='#335471')
    axs.plot(tnoev/1e6,noevData[0]['peri']*180/np.pi,'.',color=simulationcolor,rasterized=True)
    return fig,axs