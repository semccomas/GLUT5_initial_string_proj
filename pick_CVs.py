#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 15:38:20 2020

@author: semccomas
"""
import numpy as np
## so I did a dumb thing and my output array just appends the values => 0:5000 is PCA, 5000:10000 is RF...

feature_to_resids = np.load('samples_features_arr/feature_to_resids.clipped.npy')


#choose which results you want, I have done demystify with both ca distances and all heavy atom distances
atom_type = 'ca'
#atom_type = 'heavy'


def find_tops(name, cutoff):
    arr = np.load('../demystifying-master/output_files_SM/top_50_features_%s.%s.npy' %(name, atom_type))
    unique, counts = np.unique(arr[:,0], return_counts=True)
    top_counts = unique[np.where(counts >= cutoff)]
    return top_counts

pca = find_tops('PCA', 100)
rf = find_tops('RF', 20)
kl = find_tops('KL', 100)
mlp = find_tops('MLP', 90)




#%%

######################################################################
######## COMPARE TOP HITS BETWEEN METHODS TO FIND CONSISTENCIES #######

def compare_tops(name1, name2, min_counts):
    matches = []
    print('minimum # repeats saved == %i' %min_counts)
    arr1 = np.load('../demystifying-master/output_files_SM/top_50_features_%s.%s.npy' %(name1, atom_type))[:,0]
    arr2 = np.load('../demystifying-master/output_files_SM/top_50_features_%s.%s.npy' %(name2, atom_type))[:,0]
    unique1, counts1 = np.unique(arr1, return_counts=True)
    unique2, counts2 = np.unique(arr2, return_counts=True)
    for feature in unique1:
        if feature in unique2:
            c1 = counts1[np.where(unique1 == feature)[0][0]]  #np.unique sorts these features so we need to index by where the feature is
            c2 = counts2[np.where(unique2 == feature)[0][0]]
            if c1 >= min_counts and c2 >= min_counts:
                matches.append(list(feature_to_resids[int(feature)]))
                print('feature %i (resid %s) matches with %i counts in %s, %i counts in %s' %(feature,
                                                    feature_to_resids[int(feature)], c1, name1, c2, name2))
            

   
    print()
    print()
    
    return matches

min_counts = 5   ## == minimum number of repeats to compare (so like finds_tops but a much
### lower threshold because you are comparing two methods)
#a= compare_tops('PCA', 'RF', min_counts)
#compare_tops('PCA', 'KL', min_counts)
#compare_tops('PCA', 'MLP', min_counts)
#compare_tops('MLP', 'RF', min_counts)
#compare_tops('MLP', 'KL', min_counts)
#compare_tops('KL', 'RF', min_counts)






##############################################################################
############## REGION OF INTEREST IDENTIFICATION ############################
### identify residue pairs for regions of interest, using results from finding top # repeats

methods = [rf, pca, mlp, kl]
method_names = ['rf', 'pca', 'mlp', 'kl']
#important_res = np.arange(289, 296)   #tm7
#important_res = np.arange(24, 40)   #tm1
#important_res = np.arange(134, 145)   #tm4
important_res = np.arange(385, 395)   #tm10

for met, met_name in zip(methods, method_names):
    for n in important_res:
        for feat in met:
            if n in feature_to_resids[int(feat)]:   ## if this residue is found in the top # repeats
            ## but you have to relate to feature_to_resids first
                print(feature_to_resids[int(feat)], met_name, feat)




############################################################################
############### MAPPING TOP COUNTS TO RESIDUES #############################

#just take dictionary values from the keys from the 'tops'
pca_res = np.vstack([feature_to_resids[x] for x in pca.astype(int)])
kl_res = np.vstack([feature_to_resids[x] for x in kl.astype(int)])
rf_res = np.vstack([feature_to_resids[x] for x in rf.astype(int)])
mlp_res = np.vstack([feature_to_resids[x] for x in mlp.astype(int)])




