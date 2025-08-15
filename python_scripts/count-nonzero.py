import sys
import os
script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")
import matplotlib.colors as mcolors
import matplotlib
from al26_plot import read_yields,read_state,calc_cdf,use_tex
from al26_nbody import Yields,State,Metadata
from glob import glob
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.ticker as mticker
import matplotlib.cm as cm
from tqdm import tqdm
import seaborn as sns
import pandas as pd
from amuse.units import units


df = pd.read_pickle("all-sims-ratios.pkl.zst")
df = df[df.mass > 0.1]
df = df[df.mass <= 13.0]

models = ['local+sne']
isotopes = df["isotope"].unique()
nstars = [100,1000,10000]
rcs = df["rc"].unique()
nsims = sorted(df["sim_number"].unique())

print(nsims)

results = {}
results["rc"] = []
results["nstars"] = []
results["count"] = []

for isotope in isotopes:
  dfi = df[df.isotope == isotope]
  for model in models:
    dfm = dfi[dfi.model == model]
    for rc in rcs:
      dfr = dfm[dfm.rc == rc]
      for nstar in nstars:
        dfn = dfr[dfr.nstars == nstar]
        for sim in nsims:
          dfcut = dfn[dfn.sim_number == sim]
          dfcutcount = dfcut.shape[0]
          

          dfcutnonzero = dfcut[dfcut.yield_ratio_decay != 0.0]

          print(dfcut[dfcut.yield_ratio_decay != 0.2*5.25e-5])

          dfcutnonzerocount = dfcutnonzero.shape[0]
          # print(isotope,model,rc,nstar,sim)
          # print(dfcutcount,dfcutnonzerocount)
          print(dfcutnonzerocount,dfcutcount)
          if dfcutcount != 0:
            count = dfcutnonzerocount/dfcutcount

            results["rc"].append(rc)
            results["nstars"].append(nstar)
            results["count"].append(count)

df = pd.DataFrame.from_dict(results)

print(df)

use_tex()
fig = plt.figure(figsize=(7.5,3))

sns.stripplot(data=df,x="rc",y="count",hue="nstars",dodge=True,jitter=0.33,palette=sns.color_palette(),order=["0.3","1.0","3.0"],linewidth=1,alpha=0.75)

# sns.boxplot(data=df,x="rc",y="count",hue="nstars",order=["0.3","1.0","3.0","10.0"],palette=sns.color_palette())

plt.xlabel("Star-forming region radius, $R_\\textrm{c}$")
plt.ylabel("Disk enrichment fraction, $Z_\\textrm{en}$")

plt.ylim(0,1)

plt.legend(title="$N_{\\star}$")

plt.grid(ls=":")


plt.savefig("zero-yield-counts.pdf",bbox_inches="tight")


