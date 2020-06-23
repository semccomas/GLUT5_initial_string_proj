#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:26:36 2020

@author: semccomas
"""

import mdtraj as md
import numpy as np

traj = md.load_xtc('../input_f/GLUT5_occ.xtc', top = '../input_f/GLUT5_occ.start.gro')
#contacts = md.compute_contacts(traj)

np.save('samples.npy', contacts[0])
np.save('feature_to_res_ids.npy', contacts[1])