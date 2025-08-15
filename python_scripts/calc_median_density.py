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
from al26_plot import read_state,calc_cdf,use_tex,calc_local_densities_nonumba
import pandas as pd
import seaborn as sb
from math import floor,ceil
from numba import njit,prange

use_tex()
plt.figure(figsize=(5,5))

data = {}
data["nstars"]        = []
data["rc"]            = []
data["sim_number"]    = []
data["time"]          = []
data["star"]          = []
data["local_density"] = []

simsets = sorted(glob("./pt-*/pt-*-*/"))
tt1 = tqdm(simsets,position=0)
for simset_n,simset in enumerate(tt1):
  tt1.set_description(simset)
  sims = sorted(glob(simset+"pt-*/"))
  tt2 = tqdm(sims,position=1,leave=False)
  for sim_number,sim in enumerate(tt2):
    state_fnames = sorted(glob(sim+"*-state-*.zst"))
    tt3 = tqdm(state_fnames,position=2,leave=False)
    for state_fname in tt3:
      state = read_state(state_fname)
      cluster  = state.cluster
      metadata = state.metadata
      t = metadata.time.value_in(myr)
      rc = metadata.args.rc
      nstars = metadata.args.n
      local_densities = calc_local_densities_nonumba(cluster)

      # Add simulation data to dictionary
      for star_number,star in enumerate(cluster):
        density = local_densities[star_number]
        data["nstars"].append(nstars)
        data["rc"].append(rc)
        data["sim_number"].append(sim_number)
        data["time"].append(t)
        data["star"].append(star_number)
        data["local_density"].append(density)

# Write final file
df = pd.DataFrame.from_dict(data)
df.to_pickle("local-density-data.pkl.zst")
print("Finished processing!")
  
