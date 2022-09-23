from astropy.table import Table, vstack
import numpy as np
import argparse
import os
import gzip
from astropy.io import fits


tab1 = Table.read("16r-emitters/Halpha-DR3-SPLUS-PStotal-STAR-13r16-v5-Final-merge.ecsv", format="ascii.ecsv")
tab2 = Table.read("16r18-emitters/Halpha-DR3-SPLUS-PStotal-STAR-16r18-v5-Final-merge.ecsv", format="ascii.ecsv")
tab3 = Table.read("18r20-emitters/Halpha-DR3-SPLUS-PStotal-STAR-18r20-v5-v2-Final-merge.ecsv", format="ascii.ecsv")
tab4 = Table.read("20r21-emitters/Halpha-DR3-SPLUS-PStotal-STAR-20r21-v5-Final-merge.ecsv", format="ascii.ecsv")


# Merge the tables
table_merge = vstack([tab1, tab2, tab3, tab4])
print("#############################################")    
print("Number of total of objects:", len(table_merge))

table_merge.write("Halpha-DR3_PStotal-STAR_total-clean.ecsv", format="ascii.ecsv", overwrite=True)

# Save dataframe
df = (table_merge.to_pandas())
df.to_csv("Halpha-DR3_PStotal-STAR_clean.csv", index=False)
