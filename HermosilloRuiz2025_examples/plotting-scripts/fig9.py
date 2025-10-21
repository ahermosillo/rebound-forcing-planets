import matplotlib.pyplot as plt
import numpy as np
import glob
plt.style.use('paper_labels_colors.mplstyle')
from plotting import *



taus = np.array(([0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0],
[0.0,5.13e6,0.0,0.0,0.0],
[1.2e7, 5.13e6, 0, 0, 0]))*2*np.pi

fparams = np.array(([0.0,0.0,0.0,26.47],
                    [0.0,0.0,0.04,0.022],
                    np.array([0.0,0.0,0.0,0.0])*np.pi/180,
                    np.array([0.0,0.0,0.0,0.0])*np.pi/180,
                    np.array([0.0,0.0,0.0,0.0])*np.pi/180))

expEq = np.array(([None,None,None,None,None],
[None,None,None,None,None],
[None,0,None,None,None],
[0,0,None,None,None]))

indxPlanets = [2,3]
file = "../4planet-Urn-Nep-edamp-fig9/fourplanet-urn-nep-dampboth.txt"
dat = get_data(file,4,rebound=True)
fig,axs=publication_plots_ae(dat,indxPlanets,fparams,taus,expEq,rebound=True,verbose=False)
fig.supxlabel("time (Myr)")
fig.subplots_adjust(hspace=0.1,top=0.98,right=0.98,left = 0.15,bottom=0.17,wspace=0.4)
fig.savefig(f"fig9.png",dpi=300)
plt.close()