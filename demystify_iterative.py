#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:30:20 2020

@author: semccomas
"""
## taken from demo.py in demystifying

from __future__ import absolute_import, division, print_function

import logging
import numpy as np
from demystifying import feature_extraction as fe, visualization
from demystifying.data_generation import DataGenerator

logger = logging.getLogger("demo")
logger.setLevel('INFO')

num_runs = 10
outf = 'RF_%i_runs' %num_runs

samples = np.load('../scripts/samples.ca.npy')
feature_to_resids = np.load('../scripts/feature_to_res_ids.ca.npy')
labels = np.load('../scripts/labels_all.npy')

for run in range(num_runs):
    extractor = fe.RandomForestFeatureExtractor(samples=samples, labels=labels)
    extractor.extract_features()
    
    postprocessor = extractor.postprocessing()
    postprocessor.average()
    postprocessor.persist()
    
        