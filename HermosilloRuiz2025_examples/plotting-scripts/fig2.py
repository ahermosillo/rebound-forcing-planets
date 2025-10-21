import matplotlib.pyplot as plt
import numpy as np
import glob
plt.style.use('paper_labels_colors.mplstyle')
from plotting import *


### two planet example: now instead of evolving ALL ELEMENTS decided to just do exponential for a, e, i.
taus = np.array(([1e7, 5e6, 2e7, 0.0, 0.0],[1e7, 5e6, 2e7, 0.0, 0.0]))*2*np.pi
fparams = np.array(([5.0,30.0],  
                    [0.1,0.3],
                    np.array([2.0,2.0])*np.pi/180,
                    np.array([0.0,0.0])*np.pi/180,
                    np.array([0.0,0.0])*np.pi/180))
expEq = np.array(([0,0,0,None,None],[0,0,0,None,None]))

indxPlanets = [0,1]
file = "../2planet-aei-edamp-fig2/twoplanet-edamp-1e9yrs.txt"
dat = get_data(file,2,rebound=True)
fig,axs = publication_plots(dat,indxPlanets,fparams,taus,expEq,rebound=True,verbose=False) # this just plots a, e, i 
axs[0].set_xlim(0,150)
fig.supxlabel("time (Myr)")
fig.subplots_adjust(hspace=0.1,top=0.98,right=0.98,left = 0.15,bottom=0.17,wspace=0.4)
fig.savefig(f"fig2.png",dpi=300)
plt.close()