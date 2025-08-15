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
isotopes = ["26al"]
nstars = [100,1000,10000]
rcs = [0.3,1.0,3.0]
nsims = sorted(df["sim_number"].unique())

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

        ratios = []
        # count_sol = 0

        for sim in nsims:
          dfcut = dfn[dfn.sim_number == sim]
          dfcut = dfcut[dfcut.initial_mass <= 3.0]
          dfcut = dfcut[dfcut.initial_mass >= 0.1]
          dfcutcount = dfcut.shape[0]

          dfsolar = dfcut[dfcut.yield_ratio_decay >= 0.1*5.85e-5]
          dfsolarcount = dfsolar.shape[0]

          # count_sol = 0
          if dfcutcount != 0:
            # count_tot += dfcutcount
            # count_sol += dfsolarcount
            ratios.append(dfsolarcount/dfcutcount)

          

          
        
        print(rc,nstar)
        ratios = np.asarray(ratios)
        mn = ratios.mean()
        err = ratios.std() * (len(ratios))**-0.5
        print("${:.3f} \\pm {:.3f}$".format(mn,err))
        # print(count_sol,count_tot,count_sol / count_tot)