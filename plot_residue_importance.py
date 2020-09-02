#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 16:23:07 2020

@author: semccomas
"""

import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md


### transform to actual residues from clipped
pdb = md.load_pdb('../input_f/protein_only/GLUT5_in.clipped.protein.start.pdb')
top = pdb.topology
resid_dict = {}
for resid in top.residues:
    resid_dict[resid.index] = resid.resSeq

## comparing methods
def transform_importance(name):
    new = np.zeros(resid_dict[len(resid_dict)-1])   ## should be 455, is as long as # residues
    arr = np.load('../demystifying-master/closest_ca/%s/importance_per_residue.npy' %name)
    for index, importance in enumerate(arr):
        new[resid_dict[index]] = importance
    return new

pca = transform_importance('PCA')
rf = transform_importance('RF')
kl = transform_importance('KL')
mlp = transform_importance('MLP')

plt.plot(pca, color = '#0958E9', label = 'PCA')
plt.plot(rf, color = '#E99A09', label = 'RF')
plt.plot(kl, color = '#09E997', label = 'KL')
plt.plot(mlp, color = '#B5CFFF', label = 'MLP')

plt.legend()
plt.xlim(0)
plt.ylim(0)
plt.show()

plt.clf()



## comparing states. Printing residue that you find importance and then you can use methods from pick_CV's to find out
## what these residues are matching to 
rf = np.load('../demystifying-master/closest_ca/RF/importance_per_residue_and_cluster.npy')
kl = np.load('../demystifying-master/closest_ca/KL/importance_per_residue_and_cluster.npy')
mlp = np.load('../demystifying-master/closest_ca/MLP/importance_per_residue_and_cluster.npy')

cutoff = 0.6
for state, importance in enumerate(mlp.T):
    for residue in np.where(importance > cutoff)[0].tolist():
        print(resid_dict[residue], state)  
        
def transform_importance_plot_per_state(name):
    new = np.zeros((resid_dict[len(resid_dict)-1],5))   ## should be 455, is as long as # residues
    arr = np.load('../demystifying-master/closest_ca/%s/importance_per_residue_and_cluster.npy' %name)
    for state, values in enumerate(arr.T):
        for index,importance in enumerate(values):
            new[int(resid_dict[index]), state] = importance    # at this column (state #), change from zero to importance at that residue
    
    plt.plot(new[:,0], label = 'in', color = '#DBE0E9')
    plt.plot(new[:,1], label = 'in occ', color = '#95B2E9')
    plt.plot(new[:,2], label = 'occ', color = '#4F85E9')
    plt.plot(new[:,3], label = 'out occ', color = '#0958E9')
    plt.plot(new[:,4], label = 'out', color = '#000C60')
    
    plt.legend()
    plt.title(name)
    plt.xlim(0,455)
    plt.ylim(0)
    plt.xlabel('Residue')
    plt.ylabel('Importance')
    plt.show()
    return new

a = transform_importance_plot_per_state('KL')



