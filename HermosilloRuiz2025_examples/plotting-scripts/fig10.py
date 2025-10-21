import matplotlib.pyplot as plt
import numpy as np
import glob
plt.style.use('paper_labels_colors.mplstyle')
from plotting import *

dat1 = get_data("../1planet-with-without-vusr-fig10/oneplanet-edamp-novusr.txt",1,rebound=True)
dat2 = get_data("../1planet-with-without-vusr-fig10/oneplanet-edamp-yesvusr.txt",1,rebound=True)

#### one planet edamp no vusr.  change this stuff. 
taus = np.array(([1e7, 1e7, 1e7, 1e7,1e7],[1,2,3,4,5]))*2*np.pi
fparams = np.array(([7.0],
                    [0.1],
                    np.array([15.0])*np.pi/180,
                    np.array([85.0])*np.pi/180,
                    np.array([55.0])*np.pi/180))
expEq = np.array(([None,0,None,None,None],[0,0,0,0,0]))


fig,axs = plot_vusr_novusr_ex(dat1,dat2,fparams,taus,expEq)

axs[0].set_ylim(4.6,5.1)
axs[1].set_ylim(-0.05,0.35)
axs[2].set_ylim(9,11)
axs[4].set_ylim(25,35)
fig.supxlabel("time (Myr)")
fig.subplots_adjust(hspace=0.1,top=0.98,right=0.98,left = 0.15,bottom=0.17,wspace=0.4)
fig.savefig("fig10.png", dpi=300, bbox_inches='tight', metadata={"Creator": "Matplotlib"})
plt.close()