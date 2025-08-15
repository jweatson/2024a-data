import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
from glob import glob
import sys
script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")

from amuse.units import units
from al26_nbody import State,Metadata,Yields,myr,pc,msol
from al26_plot import read_state,calc_cdf,use_tex
import pandas as pd
import seaborn as sb
from math import floor,ceil

d = pd.read_pickle("all-sims-ratios.pkl.zst")
# Remove duplicates
d = d[d.isotope == "26al"]
d = d[d.model == "local"]
# Make massive star subset
df = d[d.initial_mass >= .01]

print(df)

masses = np.log10(df.initial_mass)
nbins = ceil(max(masses)) - floor(min(masses))
print(nbins)
use_tex()

plt.figure(figsize=(8,4))
# sb.histplot(masses,bins=nbins*100)
# sb.displot(masses,kind="ecdf",height=4,aspect=2.5)
sb.histplot(masses,bins=100)
plt.ylabel("Count")
plt.xlabel("Log initial mass, $\\log_{10}{(M_i)}$ ($\\log_{10}$M$_\\odot$)")
plt.title("Total star count = {}".format(len(masses)))
plt.ylim(0,50000)
plt.xlim(-2,2)
plt.grid(ls=":")

ax2 = plt.twinx()
sb.ecdfplot(masses,ax=ax2,color="tab:orange")
plt.ylabel("CDF")

# plt.yscale("log")
# plt.xscale("log")
plt.savefig("star-stats.pdf",bbox_inches="tight")

# print(d.mass)