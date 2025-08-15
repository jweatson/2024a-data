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




radii = [0.3,1.0]
nstars = [100,1000,10000]
nstars_col = ["tab:blue","tab:orange","tab:green"]
radii_style = ["-","--"]
models = ["local+sne","global+sne","sne"]


isotopes = ["26al","60fe"]
ssenrich = [[5.83e-5],[1.150e-8,1e-6]]
ssenrichl = [["$Z_{26\\mathrm{Al},\\odot} = 5.83\\times 10^{-5}$"],
             ["$Z_{60\\mathrm{Fe},\\odot} = 1.15\\times 10^{-8}$","$Z_{60\\mathrm{Fe},\\odot} = 1.00\\times 10^{-6}$"]]

for x,isotope in enumerate(isotopes):
  d = pd.read_pickle("global-model-virial-all-sim-ratios.pkl.zst")
  d = d[d.isotope == isotope]
  for model in models:
    df = d[d.model == model]
    use_tex()
    fig = plt.figure(figsize=(7.5,3))
    for i,nstar in enumerate(nstars):
      for j,radius in enumerate(radii):
        dfn = df[df.nstars == nstar]
        dfr = dfn[dfn.rc == radius]
        nstarnomcut = len(dfr)
        dfr = dfr[dfr.mass >= 0.1]
        dfr = dfr[dfr.mass < 3.0]
        nstarmcut = len(dfr)
        xx,yy = calc_cdf(dfr.yield_ratio_decay)
        label = "$r_c$={:.1f}, $N_\\star$=$10^{:.0f}$".format(radius,np.log10(nstar))
        plt.plot(xx,yy,label=label,color=nstars_col[i],ls=radii_style[j])

    plt.xscale("log")
    plt.xlabel("$^{26}$Al enrichment, $Z_{26\\mathrm{Al}}$")
    plt.ylabel("CDF")
    plt.xlim(1e-12,1e0)
    plt.ylim(0.0,1)
    plt.grid(ls=":")
    legend = fig.legend(loc="upper center",ncol=3,numpoints=1,prop={'size': 6})
    legend.get_frame().set_alpha(0)
    legend.get_frame().set_facecolor((0, 0, 1, 0.1))


    for ssn in ssenrich[x]:
      plt.axvline(ssn,c="red",ls=":")
      # plt.text(ssn, 0.01,"$Z_{26\\mathrm{Al},\\odot} = 5.83\\times 10^{-5}$",rotation=90,verticalalignment="bottom",horizontalalignment="right")
    plt.savefig("cdf-all-models-{}-{}-postproc.pdf".format(isotope,model),bbox_inches="tight")