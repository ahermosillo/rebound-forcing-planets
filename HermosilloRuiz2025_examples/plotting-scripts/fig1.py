import matplotlib.pyplot as plt
import numpy as np
import glob
plt.style.use('paper_labels_colors.mplstyle')
from plotting import *


taus = np.array(([1e7, 5e6, 4e6, 8e7, 2e7],[1,2,3,4,5]))*2*np.pi
fparams = np.array(([7.0],                          # a
                    [0.1],                          # e
                    np.array([15.0])*np.pi/180,     # i 
                    np.array([85.0])*np.pi/180,     # peri
                    np.array([60.0])*np.pi/180))    # node

expEq = np.array(([3,0,2,1,2],[0,0,2,2,2])) # functional form of evolution

indxPlanets = [0]
file = "../1planet-wacky-fig1/oneplanet-wacky.txt"
dat = get_data(file,1,rebound=True)
fig, axs = plot_planets(dat,fparams,taus,expEq,force_trange=None,rebound=True,verbose=False)
fig.supxlabel("time (Myr)")
fig.subplots_adjust(hspace=0.1,top=0.98,right=0.98,left = 0.15,bottom=0.17,wspace=0.4)
fig.savefig(f"fig1.png",dpi=300)
plt.close()