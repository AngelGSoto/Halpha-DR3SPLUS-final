'''
This script creates very nice jointplot (sns).
Author: Luis A. Gutiérrez-Soto
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.table import Table, vstack
import hdbscan
from pathlib import Path
import seaborn as sns
import glob
import json
from matplotlib.colors import PowerNorm
sns.set_color_codes()

# Read the astropy tables
table = Table.read("Halpha-DR3_PStotal-STAR_total-clean-Final-hdbscan.ecsv", format="ascii.ecsv")

# Coverting in dataframe tables
df = table.to_pandas()
print("Number of objects for:", len(df))

# Calculatind the colors
ri = df["r_PStotal"] - df["i_PStotal"]
rj660 = df["r_PStotal"] - df["J0660_PStotal"]


# Create DataFrame with the new colums with the errors on the colours
colum_ri = pd.DataFrame(ri, columns=['r - i'])
colum_rh = pd.DataFrame(rj660, columns=['r - J0660'])
df_ = pd.concat([df, colum_ri, colum_rh], axis=1)

# Mask to high probabilites to belong
mask_blue = df_["P(Blue)"] > df_["P(Red)"]
mask_red = df_["P(Red)"] > df_["P(Blue)"]


df_blue = df_[mask_blue]
df_red = df_[mask_red]

# Blue
nb = len(df_blue)
tipe_b = np.linspace(0, 0, num=nb)

df_blue['Groups'] = np.array(tipe_b)

# Red
nr = len(df_red)
tipe_r = np.linspace(1, 1, num=nr)

df_red['Groups'] = np.array(tipe_r)

# concatenate
df_final = pd.concat([df_blue, df_red])

print("''''''''''''''''''''''''''''''''''''''''''''''''''''''''")
print("'''' Number total of object:", len(df_final), "''''''''''''''''''''''")
print("''''''''''''''''''''''''''''''''''''''''''''''''''''''''")

# Plot
sns.set(rc={'axes.labelsize':30,
            'figure.figsize':(30.0, 30.0),
            'xtick.labelsize':15,
            'ytick.labelsize':15})

sns.set_style('ticks')
# rename column to put the right name on the lejend
#df_final.rename(columns = {'Label_hier':'Groups'}, inplace = True)


labels = ["Blue", "Red"]
g = sns.jointplot(data=df_final.replace({"Groups": {i: label for i, label in enumerate(labels)}}),
                  x="r - i", y="r - J0660", hue="Groups", height=5, ratio=2, space=0, 
                  kind="kde", shade=False, thresh=0.05, 
                   alpha=.5, marginal_ticks=True, palette={'#d7191c', '#2b83ba'}, legend = True,
                 ) # s=80, linewidth=1,


g.ax_marg_x.set_xlim(-1.5, 2)
g.ax_marg_y.set_ylim(-1, 4)

g.set_axis_labels(r'$r - i$', r'$r - J0660$', fontsize=15)

g.plot_joint(sns.kdeplot, shade=True, thresh=0.05, alpha=.5,)
ax = plt.gca()
#ax.legend(["Blue", "Red"])
g.plot_marginals(sns.kdeplot, color='b', shade=True, alpha=.2, legend=False)

g.ax_marg_x.set_facecolor('#ccffccaa')
g.ax_marg_y.set_facecolor('#ccffccaa')

plt.tight_layout()
plt.savefig("../../paper/Figs3/class-ri-rj0660-hdbscan.pdf")
