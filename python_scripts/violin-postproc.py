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


d = pd.read_pickle("global-model-virial-all-sim-ratios.pkl.zst")
# Filter out low and high mass stars
d = d[d.mass >= 0.1]
d = d[d.mass <= 3.0]

decaysen = ["disk progression","no disk progression"]
modelsen = ["SNE model","Local model \& SNe", "Local model","Global model \& SNe","Global model"]

for ii,decay in enumerate(["decay","nodecay"]):
  df = d
  df["log_yield"] = np.log10(df["yield_ratio_{}".format(decay)])
  # Remove nans
  df = df.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
  df = df[df.nstars != 1500]
  for jj,model in enumerate(["sne","local+sne","local","global+sne","global"]):
    use_tex()
    fig,ax = plt.subplots(ncols=1,nrows=2,sharex=True,figsize=(7.5,5))
    df2 = df[df.model == model]
    dfal = df2[df2.isotope == "26al"]
    dffe = df2[df2.isotope == "60fe"]

    ax[0].grid(True,ls=":",alpha=0.5)
    ax[1].grid(True,ls=":",alpha=0.5)

    sns.violinplot(ax=ax[0],data=dfal,x="rc",y="log_yield",hue="nstars",dodge=True,gap=0.1,inner="quart",palette=sns.color_palette(),cut=0,bw_adjust=.5,linewidth=1,order=["0.3","1.0","3.0"])
    sns.violinplot(ax=ax[1],data=dffe,x="rc",y="log_yield",hue="nstars",dodge=True,gap=0.1,inner="quart",palette=sns.color_palette(),cut=0,bw_adjust=.5,linewidth=1,order=["0.3","1.0","3.0"])

    ax[0].set_xlabel("")
    ax[0].set_ylabel("$^{26}$Al enrichment, $Z_{26\\mathrm{Al}}$")
    ax[0].axhline(np.log10(5.85e-5),c="red",ls=":")

    ax[1].set_xlabel("Star-forming region radius, R$_\\mathrm{c}$ (pc)")
    ax[1].set_ylabel("$^{60}$Fe enrichment, $Z_{60\\mathrm{Fe}}$")
    ax[1].axhline(np.log10(1.150e-8),c="red",ls=":")
    ax[1].axhline(np.log10(1e-6),c="red",ls=":")


    # Use common legend and log scale
    for axes in ax:
      axes.get_legend().remove()
      axes.set_ylim(-18,-1)
      axes.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
      ymin, ymax = axes.get_ylim()
      tick_range = np.arange(np.floor(ymin), ymax,2)
      axes.yaxis.set_ticks(tick_range)
      axes.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)], minor=True)
    lines_labels = [ax[0].get_legend_handles_labels()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    fig.legend(lines, labels,loc="upper center",ncol=4,title="$N_{\\star}$")
    plt.title("{}, {} simulated".format(modelsen[jj],decaysen[ii]))

    # Plot figure!
    outfilename = "violin-{}-{}-postproc-merge.pdf".format(model,decay)
    plt.savefig(outfilename,bbox_inches="tight")
    print("Finished "+outfilename)