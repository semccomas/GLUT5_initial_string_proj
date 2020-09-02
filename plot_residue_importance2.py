#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 13:43:59 2020

@author: semccomas
"""


import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import sklearn.cluster as sk

print('you can change both # bins as well as # clusters here, good to play around and see what you get!')
print('less bins = more features included')
print('more clusters = more separation')

### transform to actual residues from clipped
pdb = md.load_pdb('../input_f/protein_only/GLUT5_in.clipped.protein.start.pdb')
feature_to_resids = np.load('../scripts/samples_features_arr/feature_to_resids.clipped.ca.npy')
top = pdb.topology
resid_dict = {}
for resid in top.residues:
    resid_dict[resid.index] = resid.resSeq

def transform_importance(iteration, name):
    new = np.zeros(resid_dict[len(resid_dict)-1])   ## should be 455, is as long as # residues
    arr = np.load('../demystifying-master/iterative_RF/iteration%i/%s/importance_per_residue.npy' %(iteration, name))
    for index, importance in enumerate(arr):
        new[resid_dict[index]] = importance
    return new

rf_res = []
rf_feat = []
for iteration in range(1,100):
    rf_res.append(transform_importance(iteration, 'RF'))
    
    feature_importance = np.load('../demystifying-master/iterative_RF/iteration%i/RF/feature_importance.npy' %iteration)
    rf_feat.append(np.mean(feature_importance, axis = 1))

rf_res = np.vstack(rf_res)
rf_feat = np.vstack(rf_feat)

############################## 
### plotting   ################
###############################

## residue importance profile
m_rf_res = np.mean(rf_res, axis = 0)
s_rf_res = np.std(rf_res, axis = 0)
plt.plot(m_rf_res, color = 'green', linewidth = 1)
plt.fill_between(np.arange(len(m_rf_res)), m_rf_res+s_rf_res, m_rf_res-s_rf_res, color = 'green', alpha = 0.4)
plt.ylim(0)
plt.title('Residue importance')
plt.xlabel('Residue num')
plt.ylabel('Importance')
plt.show()


plt.clf()

## feature importance profile
plt.figure(1, figsize= (2,10))
m_rf_feat = np.mean(rf_feat, axis = 0) 
s_rf_feat = np.std(rf_feat, axis = 0)
plt.plot(m_rf_feat, linewidth = 0.2)
plt.ylim(0)
plt.xlabel('Feature num')
plt.ylabel('Importance')
plt.title('Feature importance')
plt.show()


plt.clf()


### histogramming
### for the bins, you can make them as fat or as skinny as you want, just change bin_number

m_rf_feat_logtrans = np.log10(m_rf_feat + 0.0000001)   #some importances are zero, still need these, so make really small instead
bin_number = 30
custom_bins = np.linspace(np.min(m_rf_feat_logtrans), np.max(m_rf_feat_logtrans), bin_number) ## start at lowest
## possible value, end in highest, make as many bins as bin_number
feat_logtrans_histogram = plt.hist(m_rf_feat_logtrans, bins = custom_bins)
max_features = np.where(m_rf_feat_logtrans > feat_logtrans_histogram[1][-2])[0]   #finding the index
### of the features that are at the rightmost edge of the histogram (the feat_log_histogram[-1] will be 
### the top value, so nothing can be above that, need to take the secondmost top value)
plt.show()


plt.clf()


#plotting max_features onto original features plot, to be sure we have good coverage
plt.plot(m_rf_feat, color = 'grey', linewidth = 0.2)
zeros = np.zeros(np.shape(m_rf_feat))
zeros[max_features] = m_rf_feat[max_features]  #just masking over values that aren't in max_features
plt.plot(zeros, color = 'blue', linewidth = 0.2)
plt.title('Feature importance, highlighted features')
plt.show()


plt.clf()


### show distribution of max_features as residue plot
max_res = feature_to_resids[max_features]
plt.scatter(max_res[:,0], max_res[:,1]) ## plot scatter of features, by x and y
plt.xlim(0, 455)
plt.ylim(0, 455)
plt.title('Features scattered')
plt.show()
plt.clf()


#################################################
########## clustering results from histogram ####
################################################

### first make elbow plot to determine proper # clusters 
score = []
for i in range(1, 20):
    kmeans = sk.KMeans(n_clusters = i, init='k-means++', max_iter = 300, n_init = 10, random_state = 0)
    kmeans.fit(max_res)
    score.append(kmeans.inertia_)

plt.plot(range(1,20), score)
plt.show()

plt.clf()


## actual clustering
n_clusters = 14
kmeans = sk.KMeans(n_clusters = n_clusters, init='k-means++', max_iter = 300, n_init = 10, random_state = 0)
pred_y = kmeans.fit_predict(max_res) 
colormap = 'plasma'
plt.scatter(max_res[:,0], max_res[:,1], c = kmeans.labels_, cmap = colormap)   #plot residue plot again, labels are indexed same as residues, so can color by cluster 
plt.colorbar()
plt.scatter(kmeans.cluster_centers_[:,0], kmeans.cluster_centers_[:,1], s=20, c = 'red', edgecolor = 'black')#np.arange(kmeans.n_clusters), cmap = colormap)  ## overlay cluster center onto residue plot
plt.xlim(0,470)
plt.ylim(0,470)
plt.savefig('../images_figures/cluster_plots/%i_clusters.%i_bins.png' %(n_clusters, bin_number), dpi = 500)

plt.show()


## reporting cluster groups
for res_cluster in range(0, n_clusters):    
    print('cluster number %i' %res_cluster)
    cluster_location = np.where(kmeans.labels_ == res_cluster)[0]
    for c in cluster_location:
        print('residue pair %s = feature %s. \n Importance %.4f and std %.4f' %(max_res[c], max_features[c], m_rf_feat[max_features[c]], s_rf_feat[max_features[c]]))  #labels is as long as max_clusters, so you can just index
        print()
        ### max_res and max_features are the same variables in two different descriptions, features is the feature index, res is the residues corresponding to feature
        ### index. So you can then use max_features index to index again onto m_rf_feat (so #18 in max_features is == features 26811. Since m_rf_feat is an array of avg
        ### importances, you can index with 26811 and find corresponding importance)
    print()
    print()
    print()

print('cluster centers = %s' %kmeans.cluster_centers_)





'''
#### old stuff 
high_values = np.where(m_rf_feat > 0.045)
important_clusters = [(0,5509), (13680,18600), (26216,27763), (32762,35777), (36768,36900), (39506,42243)]
index = 0 

for index in range(len((important_clusters))):
    print('Detecting for group %s, == approx resids %s to %s' % (important_clusters[index], feature_to_resids[important_clusters[index][0]][0],
                                                                                feature_to_resids[important_clusters[index][1]][0]))
    pair2 = []
    for v in high_values[0]:
        if v > important_clusters[index][0] and v < important_clusters[index][1]:
            print('Feature %s == resid group %s' %(v, feature_to_resids[v]))
            pair2.append(feature_to_resids[v][1]) #this is usually repeating, want to find repeats
    pair2 = np.array(pair2)
    unique, counts = np.unique(pair2, return_counts=True)
    print('Common pairing residues are %s, with %s counts' %(unique, counts))
    print()
    print()
    print()

    
'''