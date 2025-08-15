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
from al26_plot import read_state,calc_cdf,use_tex,calc_local_densities
import pandas as pd
import seaborn as sb
from math import floor,ceil
from numba import njit,prange

d = pd.read_pickle("local-density-data.pkl.zst")

nstars = [1000]
rcs = [0.3]

use_tex()
plt.figure(figsize=(5,3))

for nstar in nstars:
  for rc in rcs:
    for sim in range(32):
      print(nstar,rc)
      df = d
      df = df[df.sim_number == sim]
      df = df[df.rc == rc]
      df = df[df.nstars == nstar]
      # print(df)
      # sys.exit()
      
      times = sorted(df.time.unique())
      median_densities = []

      for t in times:
        dft = df[df.time == t]
        local_rhos = dft.local_density
        median_rho = local_rhos.median()
        median_densities.append(median_rho)

      plt.plot(times,median_densities)


plt.yscale("log")
plt.xlabel("Time, $t$ (Myr)")
plt.ylabel("Median Density, $\\rho_{1/2}$ (M$_\\odot\,$pc$^{-3}$)")
plt.savefig("median_density.pdf",bbox_inches="tight")
