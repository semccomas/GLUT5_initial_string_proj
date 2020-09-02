#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 14:46:30 2020

@author: semccomas
"""


## here I will plot the distances of the CV's as I have chosen them. I just manually write this list of chosen residue distances, taken from the output
## of plot_residue_importance2.py
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

chosen_distances = [40345, 15054, 40407, 35717, 16538, 26431, 40428, 33836, 27115, 26811, 34812, 40377, 26771, 36811]
feature_to_resids = np.load('./samples_features_arr/feature_to_resids.clipped.ca.npy')
distances = np.load('./samples_features_arr/samples.clipped.ca.normalized.npy')

a= distances[:,40345]

a = gaussian_filter(a, sigma = 2, mode = 'mirror')

plt.plot(a)