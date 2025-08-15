
import pickle
import ubjson
import zstandard as zstd
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import numba as nb
from glob import glob
import matplotlib.cm as cm
import sys,os

script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")

from al26_nbody import State,Metadata,Yields,myr,pc,msol
import scipy
from amuse.units import units

from al26_plot import calc_disk_final_enrichment,read_yields,read_state,use_tex,calc_sn_times

import pandas as pd

sim = sys.argv[1]
print(sim)

yields_fname = sorted(glob(sim+"/*yields*ubj.zst"))[-1]
state_fnames = sorted(glob(sim+"/*-state-*.zst"))
first_state_fname = state_fnames[0]
last_state_fname = state_fnames[-1]

first_state = read_state(first_state_fname)
cluster = first_state.cluster
metadata = first_state.metadata

sim_yield = read_yields(yields_fname)
lifetimes = cluster.tau_disk.value_in(myr)
sn_times,sn_masses = calc_sn_times(cluster)

print(sn_times)
print(sn_masses)
# sim_yield = calc_disk_final_enrichment(sim_yield,lifetimes)

nstars = len(cluster.mass)
cmaps = np.linspace(0,1,nstars)

use_tex()
plt.figure(figsize=(5,3))

# Using loop, much slower but slightly finer grain control
for i,star in enumerate(tqdm(cluster)):
  color = cm.get_cmap("viridis")(cmaps[i])
  stab = star.mass_27al.value_in(msol)

  # Sort lines
  t = sim_yield.time
  yl = sim_yield.local_26al[:,i]
  ys = sim_yield.sne_26al[:,i]
  # Sort pointers
  tau = lifetimes[i]
  tym = (sim_yield.local_26al_final[i] + sim_yield.sne_26al_final[i]) / stab

  yc = (yl+ys) / stab

  tcm = np.argmax(t >= tau)

  plt.plot(t[:tcm],yc[:tcm],c=color)
  plt.scatter(t[tcm-1],yc[tcm-1],c=color,marker="x")
  # plt.plot(t,yc,alpha=0.1,c=color)

  # plt.scatter(tau,tym,c=color)

for i,tt in enumerate(sn_times):
  if sn_masses[i] >= 25.0:
    c = "grey"
  else:
    c = "red"
  plt.axvline(tt,c=c,linestyle="--")

plt.yscale("log")
plt.xlim(0,20)
plt.xlabel("Time, $t$ (Myr)")
plt.ylabel("Al$_{26}$/Al$_{27}$ ratio, $Z_\\mathrm{26Al}$")
plt.grid(which='both', color='k', alpha=0.1, linestyle=':')
plt.title(metadata.filename)
plt.savefig("time-yield.pdf",bbox_inches="tight")