import glob
import json
import matplotlib.pyplot as plt
import pandas as pd
#import StringIO
from astropy.table import Table, hstack
import seaborn as sns
import argparse
import sys
import shutil

parser = argparse.ArgumentParser(
    description="""Make a table from the S-PLUS catalogs """)

parser.add_argument("source", type=str,
                    default=" teste-program",
                    help="Name of catalog, taken the prefix ")


cmd_args = parser.parse_args()
file_ = cmd_args.source + ".dat"

tab = Table.read(file_, format="ascii")

pattern1 = "SDSS-spectra/*.pdf"
file_list1 = glob.glob(pattern1)

sp = []
print(len(tab))
for i in tab:
    for j in file_list1:
        if i["ID"].split(".")[-1].replace(".", "-") == j.split("-")[-1].split(".pd")[0]:
            sp.append(j)

for objects in sp:
    shutil.copy2(objects, "SDSS-spectra-paper")
