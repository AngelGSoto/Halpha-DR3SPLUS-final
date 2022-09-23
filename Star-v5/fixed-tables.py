from astropy.table import Table
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description="""Get colored image and cut image in the r-band""")

parser.add_argument("table", type=str,
                    default="known-PN-jplus-idr",
                    help="Name of table, taken the prefix ")

parser.add_argument("--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

args = parser.parse_args()
file_ = args.table + ".ecsv"

try:
    data = Table.read(file_, format="ascii.ecsv")
except FileNotFoundError:
    file_ = args.source + ".dat"
    data = Table.read(file_, format="ascii")

id_ = []
ra = []
dec = []
for i in data:
    if i["ID"].endswith(" '"):
        id_.append(i["ID"].split("b'")[-1].split(" ")[0])
    else:
        id_.append(i["ID"].split("b'")[-1].split("'")[0])

    ra.append(i["RA"])
    dec.append(i["DEC"])

t = Table([id_, ra, dec], names=('ID', 'RA', 'DEC'))

#convert in pandas
df = (t.to_pandas())

#save
dffile= file_.replace(".ecsv", "-rightid.csv")
df.to_csv(dffile, index = False)
