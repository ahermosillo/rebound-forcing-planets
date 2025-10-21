import matplotlib.pyplot as plt
import numpy as np
import glob
plt.style.use('paper_labels_colors.mplstyle')
from plotting import *

dat1 = get_data("../2planet-omega-ev-fig4/twoplanet-Omdamp-noev.txt",2,rebound=True)
dat2 = get_data("../2planet-omega-ev-fig4/twoplanet-wdampJ2.txt",2,rebound=True)

## two planet example now only evolving the omegas
taus = np.array(([0.0, 0.0, 0.0, 3e7, 0.0],[0.0,0.0,0.0, 1e6,0.0]))*2*np.pi
fparams = np.array(([0.0,0.0],
                    [0.0,0.0],
                    np.array([0.0,0.0])*np.pi/180,
                    np.array([150.0,150.0])*np.pi/180,
                    np.array([0.0,0.0])*np.pi/180))
expEq = np.array(([None,None,None,1,None],[None,None,None,None,None]))

fig,ax=plot_w_damp(dat1,dat2,[0],fparams,taus,expEq)
ax.set_ylabel("$\omega$ (deg)")
ax.set_xlabel("time (Myr)")
ax.set_xlim(0,50)
fig.subplots_adjust(hspace=0.1,top=0.98,right=0.98,left = 0.16,bottom=0.17,wspace=0.4) 
fig.savefig("fig4.png",dpi=300)