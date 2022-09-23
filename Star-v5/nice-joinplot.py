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
table_blue = Table.read("Blue0-Halpha-DR3_PStotal-STAR_total-clean-unique.ecsv", format="ascii.ecsv")
table_red = Table.read("Red1-Halpha-DR3_PStotal-STAR_total-clean-unique.ecsv", format="ascii.ecsv")

# Coverting in dataframe tables
df_blue = table_blue.to_pandas()
df_red = table_red.to_pandas()
print("Number of objects for blue:", len(df_blue), "and red:", len(df_red))

# concatenate
df_final = pd.concat([df_blue, df_red])

# Calculatind the colors
ri = df_final["r_PStotal"] - df_final["i_PStotal"]
rj660 = df_final["r_PStotal"] - df_final["J0660_PStotal"]

# df_final["r - i"] = ri.values
# df_final["r - J0660"] = rj660.values

# print(df_final["r - J0660"])

# Create DataFrame with the new colums with the errors on the colours
colum_ri = pd.DataFrame(ri, columns=['r - i'])
colum_rh = pd.DataFrame(rj660, columns=['r - J0660'])
df_final_ = pd.concat([df_final, colum_ri, colum_rh], axis=1)


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
df_final_.rename(columns = {'Label_hier':'Groups'}, inplace = True)


labels = ["Blue", "Red"]
g = sns.jointplot(data=df_final_.replace({"Groups": {i: label for i, label in enumerate(labels)}}),
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
plt.savefig("../../paper/Figs3/class-ri-rj0660.pdf")
