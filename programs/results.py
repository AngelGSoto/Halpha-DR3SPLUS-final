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
sns.set_color_codes()
ROOT_PATH = Path("..")

# Reading the json files with synthectic photometry of the star and giant library Pickles, A. J. (1998)
def filter_mag(e, s, f1, f2, f3):
    '''
    Calculate the colors using any of set of filters
    '''
    col, col0 = [], []
    if data['id'].endswith(e):
        if data['id'].startswith(str(s)):
            filter1 = data[f1]
            filter2 = data[f2]
            filter3 = data[f3]
            diff = filter1 - filter2
            diff0 = filter1 - filter3
            col.append(diff)
            col0.append(diff0)
    
    return col, col0

def plot_mag(f1, f2, f3):
    x, y = filter_mag("Star", "", f1, f2, f3)
    for a, b in zip(x, y):
        A1.append(a)
        B1.append(b)

# Read the file
parser = argparse.ArgumentParser(
    description="""Make a table from the S-PLUS catalogs""")

parser.add_argument("fileName", type=str,
                    default="teste-program",
                    help="Name of table, taken the prefix")

cmd_args = parser.parse_args()
file_ = cmd_args.fileName + ".ecsv"

table = Table.read(file_, format="ascii.ecsv")

ra = table["RA"]
dec = table["DEC"]

icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
gal = icrs.galactic  

l_rad = gal.l.radian
l_rad[l_rad > np.pi] -= 2. * np.pi
b_rad = gal.b.radian

#latitude in degrees
b_deg = b_rad * (180/np.pi)#gal.b.degree#b_rad * (180/np.pi)

# Sintectin MS track
A1, B1 = [], []

pattern = "../../../MS_stars/*.json"
file_list = glob.glob(pattern)

for file_name in file_list:
    with open(file_name) as f:
        data = json.load(f)
        plot_mag("F0626_rSDSS", "F0660", "F0769_iSDSS")
###############################################
# New colors    ###############################
###############################################
m = (table["e_g_PStotal"] <= 0.2) & (table["e_z_PStotal"] <= 0.2)
m1 = (table["e_u_PStotal"] <= 0.2) & (table["e_g_PStotal"] <= 0.2) & (table["e_r_PStotal"] <= 0.2) 
zg = table['z_PStotal'] - table['g_PStotal']
gr = table['g_PStotal'] - table['r_PStotal']
ug = table['u_PStotal'] - table['g_PStotal']
ri = table['r_PStotal'] - table['i_PStotal']
rz = table['r_PStotal'] - table['z_PStotal']
gz = table['g_PStotal'] - table['z_PStotal']
rj0660 = table['r_PStotal'] - table['J0660_PStotal']
zj0660 = table['z_PStotal'] - table['J0660_PStotal']

##############################################
# Plots
color_map = plt.cm.Spectral_r
color_palette = sns.color_palette('Paired', 55)
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize=(15, 11))
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    plt.xlabel(r"$r - i$", fontsize=35)
    plt.ylabel(r"$r - J0660$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    maskfw = table["FWHM"] < 1000
    scat = ax.scatter(ri[maskfw], rj0660[maskfw], s=15*table["FWHM"][maskfw], edgecolor='black',
                             c=table["r_PStotal"][maskfw], alpha=0.7, zorder = 2, cmap='hot')#RdBU_r
    #pal = sns.dark_palette("magma", as_cmap=True)
    #pal = sns.cubehelix_palette(as_cmap=True)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    axx = sns.kdeplot(B1, A1, zorder = 3, cmap=pal);
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    ax.set(
      xlim=[-1.5, 2.5],
      ylim=[-1.2, 5.5])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    cb = fig.colorbar(scat, extend='both', ax=ax)#
    cb.set_label("$r[mag]$", fontsize=35)
    cb.ax.tick_params(labelsize=30)
    #Symbol size indicates outer shell radius
    plt.text(0.02, 0.95, 'Symbol size indicates FWHM',
             transform=ax.transAxes, fontsize=22)

    # main sequence and giant stars loci
    x1, y1 = 0.3, 0.3
    el = mpatches.Ellipse((x1, y1), 0.3, 0.4, angle=30, alpha=0.3)
    ax.annotate("Contour indicates main-sequence and giant stars loci",
                xy=(0.1, -0.1), xytext=(-1.3, -1.), color='black', size=22,
                zorder= 111, arrowprops=dict(arrowstyle="fancy",
                            color="0.5",
                            patchB=el,
                            shrinkB=5,
                            connectionstyle="arc3,rad=0.3",
                                                        ))
    
    plt.savefig("../../paper/Figs3/final-emitters.pdf")
    ##########################################################
    # (g - r) vs (z - g)
    fig, ax1 = plt.subplots(figsize=(15, 11))
    ax1.spines["top"].set_visible(False)  
    ax1.spines["right"].set_visible(False)
    plt.xlabel(r"$z - g$", fontsize=35)
    plt.ylabel(r"$g - r$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    ax1.set(
        xlim=[-6.8, 2.5], ylim=[-3., 5.]
      )
    #scat = ax.scatter(zg, gr, s=15*table["FWHM"], edgecolor='black',
                             #c=table["R_PStotal"], alpha=0.7, zorder = 2, cmap='RdBu_r')
    # Limiting the blue and red region
    x_new = np.linspace(-15.0, 1000, 200)
    y = -1.8*x_new + 2.2

    ax1.plot(x_new, y, color='k', zorder=100, linestyle='-.')
    density_scatter(zg[m], gr[m], ax=ax1)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    #ax1.set(
      #xlim=[-15, 15],
      #ylim=[-15, 15])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    #cb = fig.colorbar(scat,extend='both', ax=ax).set_label("$r-band$", fontsize=35)   
    plt.savefig("../../paper/Figs3/red-blue-colorObjects-gr.pdf")

    ##########################################################
    # (r - z) vs (g - r)
    fig, ax1 = plt.subplots(figsize=(12, 12))
    ax1.spines["top"].set_visible(False)  
    ax1.spines["right"].set_visible(False)
    plt.xlabel(r"$g - r$", fontsize=35)
    plt.ylabel(r"$r - z$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    ax1.set(
        xlim=[-2.1, 3.7], ylim=[-2.2, 4.]
      )
    #scat = ax.scatter(zg, gr, s=15*table["FWHM"], edgecolor='black',
                             #c=table["R_PStotal"], alpha=0.7, zorder = 2, cmap='RdBu_r')
    # Limiting the blue and red region
    x_new_ = np.linspace(-15.0, 1000, 200)
    y_ = -1.2*x_new + 1.6 

    ax1.plot(x_new_, y_, color='k', zorder=100, linestyle='-.')
    density_scatter(gr[m], rz[m], ax=ax1)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)

    #create new axes on the right and on the top of the current axes
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    # Set aspect of the main axes.
    #ax1.set_aspect(1.)

    divider = make_axes_locatable(ax1)
    # below height and pad are in inches
    ax_histx = divider.append_axes("top", 1.5, pad=0.5, sharex=ax1)
    ax_histy = divider.append_axes("right", 1.5, pad=0.5, sharey=ax1)
    
    # make some labels invisible
    ax_histx.xaxis.set_tick_params(labelbottom=False)
    ax_histy.yaxis.set_tick_params(labelleft=False)
    ax_histx.yaxis.set_tick_params(labelsize=35) 
    ax_histy.xaxis.set_tick_params(labelsize=35)

    # now determine nice limits by hand:
    binwidth = 0.05
    xymax = max(np.max(np.abs(gr[m])), np.max(np.abs(rz[m])))
    lim = (int(xymax/binwidth) + 1)*binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(gr[m], bins=bins, histtype="stepfilled", color = "violet")
    ax_histy.hist(rz[m], bins=bins, histtype="stepfilled",  color = "violet", orientation='horizontal')

    # the xaxis of ax_histx and yaxis of ax_histy are shared with ax,
    # thus there is no need to manually adjust the xlim and ylim of these
    # axis.

    ax_histx.set_yticks([0,  150])
    ax_histy.set_xticks([0,  100])
    
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    #cb = fig.colorbar(scat,extend='both', ax=ax).set_label("$r-band$", fontsize=35)   
    plt.savefig("../../paper/Figs3/red-blue-colorObjects-rz.pdf")

    ##########################################################
    # (r - j0660) vs (g - z)
    fig, ax1 = plt.subplots(figsize=(15, 11))
    ax1.spines["top"].set_visible(False)  
    ax1.spines["right"].set_visible(False)
    plt.xlabel(r"$g - z$", fontsize=35)
    plt.ylabel(r"$r - J0660$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    #ax1.set(
        #xlim=[-2, 3.8], ylim=[-3., 5.]
      #)
    #scat = ax.scatter(zg, gr, s=15*table["FWHM"], edgecolor='black',
                             #c=table["R_PStotal"], alpha=0.7, zorder = 2, cmap='RdBu_r')
    # Limiting the blue and red region
    x_new_ = np.linspace(-15.0, 1000, 200)
    y_ = -1.45*x_new_ + 1.8

    #ax1.plot(x_new_, y_, color='k', zorder=100, linestyle='-.')
    density_scatter(rz[m], rj0660[m], ax=ax1)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    #ax1.set(
      #xlim=[-15, 15],
      #ylim=[-15, 15])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    #cb = fig.colorbar(scat,extend='both', ax=ax).set_label("$r-band$", fontsize=35)   
    plt.savefig("../../paper/Figs3/red-blue-colorObjects-gz-rj0660.pdf")

    ##########################################################
    # (z - j0660) vs (z - g)
    fig, ax1 = plt.subplots(figsize=(15, 11))
    ax1.spines["top"].set_visible(False)  
    ax1.spines["right"].set_visible(False)
    plt.xlabel(r"$z - g$", fontsize=35)
    plt.ylabel(r"$z - J0660$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    #ax1.set(
        #xlim=[-2, 3.8], ylim=[-3., 5.]
      #)
    #scat = ax.scatter(zg, gr, s=15*table["FWHM"], edgecolor='black',
                             #c=table["R_PStotal"], alpha=0.7, zorder = 2, cmap='RdBu_r')
    # Limiting the blue and red region
    x_new_ = np.linspace(-15.0, 1000, 200)
    y_ = -1.45*x_new_ + 1.8

    #ax1.plot(x_new_, y_, color='k', zorder=100, linestyle='-.')
    density_scatter(zg[m], zj0660[m], ax=ax1)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    #ax1.set(
      #xlim=[-15, 15],
      #ylim=[-15, 15])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    #cb = fig.colorbar(scat,extend='both', ax=ax).set_label("$r-band$", fontsize=35)   
    plt.savefig("../../paper/Figs3/red-blue-colorObjects-zg-zj0660.pdf")

    ##########################################################
    # (g - r) vs (u - g)
    fig, ax11 = plt.subplots(figsize=(15, 11))
    ax11.spines["top"].set_visible(False)  
    ax11.spines["right"].set_visible(False)
    plt.xlabel(r"$u - g$", fontsize=35)
    plt.ylabel(r"$g - r$", fontsize=35)
    plt.tick_params(axis='x', labelsize=35) 
    plt.tick_params(axis='y', labelsize=35)
    #scat = ax.scatter(zg, gr, s=15*table["FWHM"], edgecolor='black',
                             #c=table["R_PStotal"], alpha=0.7, zorder = 2, cmap='RdBu_r')
    density_scatter(ug[m1], gr[m1], ax=ax11)
    pal = sns.cubehelix_palette(start=1, rot=0, dark=-10, light=50, reverse=True, as_cmap=True)
    #pal = sns.color_palette("Paired", 19, as_cmap=True)
    #pal = sns.color_palette("bright")
    #ax2.plot(fit_line, 0.42917 * fit_line - 0.04333, color="k", ls="--")
    #ax1.set(
      #xlim=[-15, 15],
      #ylim=[-15, 15])
    font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }
    #cb = fig.colorbar(scat,extend='both', ax=ax).set_label("$r-band$", fontsize=35)   
    plt.savefig("../../paper/Figs3/red-blue-colorObjects-ug.pdf")
###################################################################################################################
#Distribution of Halpha emitters
with sns.axes_style("ticks"):
    fig = plt.figure(figsize=(14,7))
    ax = fig.add_subplot(1,1,1, projection='aitoff')
    plt.xlabel(r'$l (Gal)$')
    plt.ylabel(r'$b (Gal)$')
    ax.xaxis.label.set_fontsize(23)
    ax.yaxis.label.set_fontsize(23)
    plt.tick_params(axis='x', labelsize=23) 
    plt.tick_params(axis='y', labelsize=23)
    #ax.scatter(l_rad, b_rad, s=1, color='black', alpha=0.2)
    density_scatter(l_rad, b_rad, ax=ax)
    ax.grid(True, linestyle='-.', linewidth=0.7)
    #plt.colorbar(image, spacing='uniform', extend='max')
    plt.savefig("../../paper/Figs3/halpha-emitters-galactic-aitoff.pdf")
    
    ##########################
    # Bar diagram
    fig1, ax1 = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
    plt.xlabel(r"$r - J0660$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    r_j0660 = [x for x in rj0660]
    g = sns.distplot(r_j0660, 
                 norm_hist=True, kde=True, ax=ax1,
                 bins=20, hist_kws=dict(range=[-3.0, 3.0], color='r')
                )
    #ax1.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    ymax = ax.get_ybound()[1]
    sns.despine()
    plt.tight_layout()
    plt.savefig("../../paper/Figs3/distribution-Halpha.pdf")
    ##########################
    # Distribution r - i color
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
    plt.xlabel(r"$r - i$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    r_i = [x for x in ri]
    sns.distplot(r_i, 
                 norm_hist=True, kde=True, ax=ax2,
                 bins=20, hist_kws=dict(range=[-3.0, 3.0], color='r')
                )
    #ax2.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    ymax = ax.get_ybound()[1]
    sns.despine()
    plt.tight_layout()
    plt.savefig("../../paper/Figs3/distribution-ri.pdf")
    #########################
    # Distribution  r-mag
    fig3, ax3 = plt.subplots(1, 1, figsize=(10, 6), sharex=True)
    plt.xlabel(r"$r$", fontsize=33)
    #plt.ylabel(r"Density", fontsize=28)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    r = [x for x in table["r_PStotal"]]
    sns.distplot(r, 
                 norm_hist=False, kde=True, ax=ax3,
                 bins=20)#, hist_kws=dict(range=[-3.0, 3.0])
                #)
    ax3.set(xlim=[14, 22])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../../paper/Figs3/distribution_r.pdf")
    #########################
    # Distribution b coordinate
    fig4, ax4 = plt.subplots(1, 1, figsize=(10, 5), sharex=True)
    plt.xlabel(r"$b(\degree)$", fontsize=33)
    plt.ylabel(r"# of sources", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    sns.distplot(b_deg, 
                 norm_hist=False, kde=False, ax=ax4,
                 bins=30,  color='purple' #hist_kws=dict(range=[-70, 70])
                 )
    #plt.axvline(x=-15.5)
    #plt.axvline(x=-41)
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../../paper/Figs3/distribution-bgalactic.pdf")
    #########################
    # r vs b (Gal)
    fig4, ax4 = plt.subplots(1, 1, figsize=(9, 12), sharex=True)
    plt.xlabel(r"$b(Gal)$", fontsize=33)
    plt.ylabel(r"r", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    density_scatter(b_rad, table["r_PStotal"], ax=ax4)
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../../paper/Figs3/bvsr.pdf")
    ###############################################################
    # Distribution z - g coordinate
    fig5, ax5 = plt.subplots(1, 1, figsize=(11, 5), sharex=True)
    plt.xlabel(r"$z - g$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    ax5.set(
      xlim=[-5.0, 2.5]
      )
    zg = [x for x in zg[m]]
    sns.distplot(zg, 
                 norm_hist=True, kde=True, ax=ax5,
                 bins=50, hist_kws=dict(range=[-6.0, 6.0],  color='r')
                )
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.tight_layout()
    plt.savefig("../../paper/Figs3/distribution-zg.pdf")
    ###############################################################
    # Distribution z - g coordinate
    fig6, ax6 = plt.subplots(1, 1, figsize=(11, 5), sharex=True)
    plt.xlabel(r"$g - r$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    ax6.set(
      xlim=[-1.5, 2.5]
      )
    gr = [x for x in gr[m]]
    sns.distplot(gr, 
                 norm_hist=True, kde=True, ax=ax6,
                 bins=80, hist_kws=dict(range=[-6.0, 6.0],  color='r')
                )
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../../paper/Figs3/distribution-gr.pdf")
     ###############################################################
    # Distribution r - z coordinate
    fig6, ax6 = plt.subplots(1, 1, figsize=(11, 5), sharex=True)
    plt.xlabel(r"$r - z$", fontsize=33)
    plt.ylabel(r"Density", fontsize=33)
    plt.tick_params(axis='x', labelsize=33) 
    plt.tick_params(axis='y', labelsize=33)
    ax6.set(
      xlim=[-1.5, 2.5]
      )
    rz = [x for x in rz[m]]
    sns.distplot(rz, 
                 norm_hist=True, kde=True, ax=ax6,
                 bins=100, hist_kws=dict(range=[-6.0, 6.0],  color='r')
                )
    #ax4.set(xlim=[-0.7, 1.8])
    #ax.legend(loc='upper left')
    sns.despine()
    plt.tight_layout()
    plt.savefig("../../paper/Figs3/distribution-rz.pdf")
    
