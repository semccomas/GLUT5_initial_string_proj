import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plotxlim = (0.5,2)
plotylim = (0.5,2)

#### plotting density - use KDE

fig, ax = plt.subplots()

def plot_gate_dist(name, col, color_map):
    ec = np.loadtxt('../gate_dists/extracellular/%s.EC.xvg' %name)[:,1]
    ic = np.loadtxt('../gate_dists/intracellular/%s.IC.xvg' %name)[:,1]
    
    sns.kdeplot(ic,ec, shade=True, shade_lowest=False, alpha=0.5,legend=True, cbar=False, ax=ax, cmap=color_map)
   # b= ax.scatter(ic,ec, color=col, label=name, alpha = 1)   #needs to be after so that the spots come on top

   
plot_gate_dist('GLUT5_out', 'green', "Greens")
plot_gate_dist('GLUT5_out_occ', 'grey', "Greys")
plot_gate_dist('GLUT5_occ', 'red', "Reds")
plot_gate_dist('GLUT5_in_occ', 'orange', "Oranges")
plot_gate_dist('GLUT5_in', 'blue', "Blues")

plt.ylabel('Extracellular gate distance (nm)')
plt.xlabel('Intracellular gate distance (nm)')


colors = ["green", "grey", "red", "orange", "blue"]
texts = ["Out Open", "Out occ.", "Occ.", "In Occ.", "In Open"]

patches = [ plt.plot([],[], marker="o", ms=10, ls="", mec=None, color=colors[i], alpha = 0.5, 
            label="{:s}".format(texts[i]) )[0]  for i in range(len(texts)) ]   #seaborn legends are weird... 
plt.legend(handles=patches)

plt.savefig('gate_distances_KDE.png', dpi =500)
#plt.show()

plt.clf()


### plotting as a function of time - use contour


import matplotlib.pyplot as plt
import numpy as np


name='GLUT5_out'
ec = np.loadtxt('../gate_dists/extracellular/%s.EC.xvg' %name)[:,1]
ic = np.loadtxt('../gate_dists/intracellular/%s.IC.xvg' %name)[:,1]
times = np.loadtxt('../gate_dists/intracellular/%s.IC.xvg' %name)[:,0] / 1000

plt.scatter(ic, ec, c=times, cmap = 'Greens')
plt.colorbar()
#plt.xlim(plotxlim)
#plt.ylim(plotylim)
plt.show()











