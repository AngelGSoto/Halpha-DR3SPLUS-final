import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import pandas as pd
import numpy as np
from astropy.table import Table
import seaborn as sns
import argparse
import sys
import os
import glob
import json
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from pathlib import Path
from density_scatter import density_scatter
import hdbscan
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.cluster.hierarchy as shc
sns.set_color_codes()
ROOT_PATH = Path("../../paper/Figs3")

# Read the file
table = Table.read("Halpha-DR3_PStotal-STAR_total-clean-unique.ecsv", format="ascii.ecsv")

# Colors
m = (table["e_g_PStotal"] <= 0.2) & (table["e_z_PStotal"] <= 0.2)
m1 =  (table["e_u_PStotal"] <= 0.2) &(table["e_g_PStotal"] <= 0.2) & (table["e_i_PStotal"] <= 0.2) 
zg = table['z_PStotal'][m] - table['g_PStotal'][m]
gr = table['g_PStotal'][m] - table['r_PStotal'][m]
ri = table['r_PStotal'][m] - table['i_PStotal'][m]
rz = table['r_PStotal'][m] - table['z_PStotal'][m]
ug = table['u_PStotal'][m1] - table['g_PStotal'][m1]
gr_ = table['g_PStotal'][m1] - table['r_PStotal'][m1]

# Create an array
X = np.array(list(zip(gr, rz)))
print("Shape:", X.shape)
# Standarized the data
X_std = StandardScaler().fit_transform(X)

# Applying HDBSCAN
clusterer = hdbscan.HDBSCAN(min_samples=15, min_cluster_size=60, prediction_data=True).fit(X_std) # 40 60
labels_h = clusterer.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels_h)) - (1 if -1 in labels_h else 0)
n_cluster0 = list(labels_h).count(0)
n_cluster1 = list(labels_h).count(1)
n_cluster2 = list(labels_h).count(2)
n_noise_ = list(labels_h).count(-1)

# Print parameters
print('##########################################################')
print('Estimated number of clusters: %d' % n_clusters_)
print('Estimated number of cluster points 0: %d' % n_cluster0)
print('Estimated number of cluster points 1: %d' % n_cluster1)
print('Estimated number of cluster points 2: %d' % n_cluster2)
print('Estimated number of noise points: %d' % n_noise_)
print('##########################################################')

# Getting the probabilities
prob = clusterer.probabilities_

# Add label to the table and making the colors
table_= table[m]

table_["Label"] = labels_h

mask0 = table_["Label"] == -1
mask1 = table_["Label"] == 0
mask2 = table_["Label"] == 1

# Making the colors
zg_0 = table_['z_PStotal'][mask0] - table_['g_PStotal'][mask0]
gr_0 = table_['g_PStotal'][mask0] - table_['r_PStotal'][mask0]
rz_0 = table_['r_PStotal'][mask0] - table_['z_PStotal'][mask0]
zg_1 = table_['z_PStotal'][mask1] - table_['g_PStotal'][mask1]
gr_1 = table_['g_PStotal'][mask1] - table_['r_PStotal'][mask1]
rz_1 = table_['r_PStotal'][mask1] - table_['z_PStotal'][mask1]
zg_2 = table_['z_PStotal'][mask2] - table_['g_PStotal'][mask2]
gr_2 = table_['g_PStotal'][mask2] - table_['r_PStotal'][mask2]
rz_2 = table_['r_PStotal'][mask2] - table_['z_PStotal'][mask2]

# Soft clustering
soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
print(soft_clusters)

table_["P(Blue)"] = soft_clusters[:,1]
table_["P(Red)"] = soft_clusters[:,0]

#Save the tables
asciifile = "Halpha-DR3_PStotal-STAR_total-clean-Final-hdbscan.ecsv" 
table_.write(asciifile, format="ascii.ecsv", overwrite=True)
df_hd = table_.to_pandas()
df_hd.to_csv("Halpha-DR3_PStotal-STAR_total-clean-Final-hdbscan.csv", index=False)

table_[mask1].write("Halpha-DR3_PStotal-STAR_total-clean-Final-hdbscan-group0.ecsv", format="ascii.ecsv", overwrite=True)
table_[mask2].write("Halpha-DR3_PStotal-STAR_total-clean-Final-hdbscan-group1.ecsv", format="ascii.ecsv", overwrite=True)

# Equation constructed from synthetic phometry
# Limiting the blue and red region
x_new = np.linspace(-15.0, 1000, 200)
y = -1.2*x_new + 1.6 

#############################################################
#Plot the results  ##########################################
#############################################################

#Build the cluster hierarchy

fig, ax = plt.subplots(figsize=(10, 7))
clusterer.condensed_tree_.plot(select_clusters=True,
                               selection_palette=sns.color_palette('RdBu', 2),
                               colorbar=True)

fig.savefig(ROOT_PATH / "cluster-hierarchy-hdbscan.pdf")
plt.clf()

##########################################################
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig, ax1 = plt.subplots(figsize=(12, 10))

ax1.fill_between(x_new, y, -100, color="k", alpha=0.1)
ax1.plot(x_new, y, c="k", zorder=11, lw=0.5)

plt.tick_params(axis='x', labelsize=38) 
plt.tick_params(axis='y', labelsize=38)

plt.xlabel(r'$g - r$', fontsize= 38)
plt.ylabel(r'$r - z$', fontsize= 38)

ax1.scatter(
        gr_0,
        rz_0,
        marker="o",
        c=sns.xkcd_rgb["grey"],
        label="Outliers",
        edgecolors="w", alpha=0.7, zorder=3
    )

ax1.scatter(
        gr_1,
        rz_1,              
        s = 70,
        alpha = 0.7,
        c=sns.xkcd_rgb["dark pink"],
        label="Red",
        edgecolors="w", zorder=4
     )

ax1.scatter(
        gr_2,
        rz_2,
        s = 70,
        alpha = 0.7,
        c=sns.xkcd_rgb["cerulean"],
        label="Blue",
        edgecolors="w",
        zorder=4
     )

sns.kdeplot(
    gr_2,
    rz_2,
    ax=ax1,
    norm=PowerNorm(0.5), zorder=6,
        cmap="Blues",
    )

sns.kdeplot(
    gr_1,
    rz_1,
    ax=ax1,
    norm=PowerNorm(0.5), zorder=6,
        cmap="Reds",
    )

ax1.legend(ncol=1, markerscale = 2,  title="Groups", title_fontsize=32, prop={'family': 'monospace', 'size': 25}, **lgd_kws)
ax1.set(xlim=[-2.1, 3.7], ylim=[-2.2, 4.])#, xscale="log", yscale="log")
ax1.text(0.25, 1.05, "HDBSCAN", fontsize=30,
                                 bbox=dict(facecolor='gray', alpha=0.2),
                                                      transform=ax.transAxes)
#ax1.set_aspect("equal")
#ax.set(xlabel=r"$z - g$", ylabel=r"$g - r$")
plt.tight_layout()
fig.savefig(ROOT_PATH / "blue-red-hdbscan.pdf")
plt.clf()

#################################
#Soft clusters   ################
#################################
soft_clusters = hdbscan.all_points_membership_vectors(clusterer)
color_palette = sns.color_palette('Paired', 12)
cluster_colors = [color_palette[np.argmax(x)]
                  for x in soft_clusters]

# Mask to high probabilites to belong
mask_red = soft_clusters[:,0] > soft_clusters[:,1]
mask_blue = soft_clusters[:,1] >= soft_clusters[:,0]

sns.set_style('ticks')
fig, ax2 = plt.subplots(figsize=(12, 10))
plt.tick_params(axis='x', labelsize=38) 
plt.tick_params(axis='y', labelsize=38)

plt.xlabel(r'$g - r$', fontsize= 38)
plt.ylabel(r'$r - z$', fontsize= 38)
ax2.set(xlim=[-2.1, 3.7], ylim=[-2.2, 4.])#, xscale="log", yscale="log")
ax2.fill_between(x_new, y, -100, color="k", alpha=0.1)
ax2.plot(x_new, y, c="k", zorder=11, lw=0.5)
#ax1.scatter(zg, gr, s=50, linewidth=0.2, c=cluster_colors, edgecolors="w", alpha=0.25)

ax2.scatter(
        gr[mask_red], rz[mask_red],
        s = 70,
        alpha = 0.7,
        c=sns.xkcd_rgb["dark pink"],
        label="Red",
        edgecolors="w", zorder=4
    )

sns.kdeplot(
    gr[mask_red], rz[mask_red],
    ax=ax2,
    norm=PowerNorm(0.5), zorder=5,
        cmap="Reds",
 )

ax2.scatter(
        gr[mask_blue], rz[mask_blue],
        s = 70,
        alpha = 0.7,
        c=sns.xkcd_rgb["cerulean"],
        label="Blue",
        edgecolors="w", zorder=6
    )

sns.kdeplot(
    gr[mask_blue], rz[mask_blue],
    ax=ax2,
    norm=PowerNorm(0.5), zorder=7,
        cmap="Blues",
 )

ax2.legend(ncol=1, markerscale = 2,  title="Groups", title_fontsize=32, prop={'family': 'monospace', 'size': 25}, **lgd_kws)
ax2.text(0.19, 1.05, "Soft Clustering for HDBSCAN", fontsize=30,
                                 bbox=dict(facecolor='gray', alpha=0.2),
                                                       transform=ax.transAxes)
#ax2.set_aspect("equal")
#ax.set(xlabel=r"$z - g$", ylabel=r"$g -a r$")
plt.tight_layout()
fig.savefig(ROOT_PATH / "blue-red-hdbscan-soft-clustering.pdf")
plt.clf()

########################################
#Dendrograms for Hierarchical Clustering

fig, ax3 = plt.subplots(figsize=(10, 7))
#plt.figure(figsize=(10, 7))
#plt.title("Customer Dendograms")
plt.xlabel('sample index',  fontsize= 25)
plt.ylabel('distance', fontsize= 25)
plt.tick_params(axis='y', labelsize=25)

# Set the colour of the cluster here:
shc.set_link_color_palette(['#d62728', '#00008B'])
with plt.rc_context({'lines.linewidth': 2.0}):
    dend = shc.dendrogram(shc.linkage(X_std, method='ward'),
                          truncate_mode='lastp',
                          p=12,  # show only the last p merged clusters
                          leaf_rotation=45.,
                          leaf_font_size=18.,
                          show_contracted=True,
                          above_threshold_color='#5F9EA0')
#dend = shc.dendrogram(shc.linkage(X, method='ward'))

# dendrogram(
#     X,
#     truncate_mode='lastp',  # show only the last p merged clusters
#     p=12,  # show only the last p merged clusters
#     leaf_rotation=90.,
#     leaf_font_size=12.,
#     show_contracted=True,  # to get a distribution impression in truncated branches
# )
plt.tight_layout()
fig.savefig(ROOT_PATH / "Customer-Dendrograms.pdf")
plt.clf()
