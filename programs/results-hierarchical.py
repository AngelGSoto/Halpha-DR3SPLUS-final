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
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.cluster.hierarchy as shc
sns.set_color_codes()
ROOT_PATH = Path("../../paper/Figs3")

table_blue = Table.read("Blue0-Halpha-DR3_PStotal-STAR_total-clean-unique.ecsv", format="ascii.ecsv")
table_red = Table.read("Red1-Halpha-DR3_PStotal-STAR_total-clean-unique.ecsv", format="ascii.ecsv")

# Making the colors
zg_blue = table_blue['z_PStotal'] - table_blue['g_PStotal']
gr_blue = table_blue['g_PStotal'] - table_blue['r_PStotal']
rz_blue = table_blue['r_PStotal'] - table_blue['z_PStotal']

zg_red = table_red['z_PStotal'] - table_red['g_PStotal']
gr_red = table_red['g_PStotal'] - table_red['r_PStotal']
rz_red = table_red['r_PStotal'] - table_red['z_PStotal']

# Equation constructed form synthetic phometry
# Limiting the blue and red region
x_new = np.linspace(-15.0, 1000, 200)
y = -1.2*x_new + 1.6 

#######################################
lgd_kws = {'frameon': True, 'fancybox': True, 'shadow': True}
sns.set_style('ticks')
fig, ax = plt.subplots(figsize=(12, 10))

ax.fill_between(x_new, y, -100, color="k", alpha=0.1)
ax.plot(x_new, y, c="k", zorder=11, lw=0.5)

plt.tick_params(axis='x', labelsize=38) 
plt.tick_params(axis='y', labelsize=38)

plt.xlabel(r'$g - r$', fontsize= 38)
plt.ylabel(r'$r - z$', fontsize= 38)

ax.scatter(
        gr_red,
        rz_red,
        s = 70,
        c=sns.xkcd_rgb["dark pink"],
        label="Red",
        edgecolors="w", alpha=0.7, zorder=10
    )

sns.kdeplot(
    gr_red,
    rz_red,
    ax=ax,
    norm=PowerNorm(0.5), zorder=11,
        cmap="Reds",
 )

ax.scatter(
        gr_blue,
        rz_blue,
        s = 70,
        c=sns.xkcd_rgb["cerulean"],
        label="Blue",
        edgecolors="w", alpha=0.7, zorder=6
    )

sns.kdeplot(
    gr_blue,
    rz_blue,
    ax=ax,
    norm=PowerNorm(0.5), zorder=7,
        cmap="Blues",
 )

ax.legend(ncol=1, markerscale = 2,  title="Groups", title_fontsize=32, prop={'family': 'monospace', 'size': 25}, **lgd_kws)
ax.set(xlim=[-2.1, 3.7], ylim=[-2.2, 4.])#, xscale="log", yscale="log")
ax.text(0.2, 0.9, "HAC", fontsize=30,
                                 bbox=dict(facecolor='gray', alpha=0.2),
                                                       transform=ax.transAxes)
#ax.set_aspect("equal")
#ax.set(xlabel=r"$z - g$", ylabel=r"$g - r$")
plt.tight_layout()

fig.savefig(ROOT_PATH / "blue-red-hierarchical.pdf")
