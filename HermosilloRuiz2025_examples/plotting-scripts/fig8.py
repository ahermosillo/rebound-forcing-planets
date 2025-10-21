import matplotlib.pyplot as plt
import numpy as np
import glob
plt.style.use('paper_labels_colors.mplstyle')
from plotting import *


taus = np.array(([0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0],
[0.0,2e7,0.0,0.0,0.0],
[1.1875004e7, 5.13081856e+06, 0.0, 0.0, 0.0]))*2*np.pi

fparams = np.array(([0.0,0.0,0.0,26.47221282],
                    [0.0,0.0,0.04,0.0225727],
                    np.array([0.0,0.0,0.0,0.0])*np.pi/180,
                    np.array([0.0,0.0,0.0,0.0])*np.pi/180,
                    np.array([0.0,0.0,0.0,0.0])*np.pi/180))
expEq = np.array(([None,None,None,None,None],
[None,None,None,None,None],
[None,0,None,None,None],
[0,0,None,None,None]))


indxPlanets = [2,3]
file = "../4planet-Nep-edamp-fig8/fourplanet-tsiganis-Nep-edamp.txt"
dat = get_data(file,4,rebound=True)
fig,axs=publication_plots_ae(dat,indxPlanets,fparams,taus,expEq,rebound=True,verbose=False)

eeq2 = get_equation(dat[2],'e',0.03-0.06,5.13081856e6*2*np.pi,0,None)
t_eq2, eeq_line = eeq2[0],np.array(eeq2[1])
axs[1].plot(t_eq2/(1e6*2*np.pi),eeq_line,color='#354E6E',linestyle='-.',rasterized=True) 

fig.supxlabel("time (Myr)")
fig.subplots_adjust(hspace=0.1,top=0.98,right=0.98,left = 0.15,bottom=0.17,wspace=0.4)
fig.savefig(f"fig8.png",dpi=300)
plt.close()