import sys
import os
script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")
import matplotlib.colors as mcolors
import matplotlib
from al26_plot import read_yields,read_state,calc_cdf,use_tex,get_high_mass_star_indices
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

# Get full list of simulation sets
simsets = sorted(glob("./pt-*/pt*/"))

tt1 = tqdm(simsets,position=0)
for simset in tt1:
  tt1.set_description(simset)
  sims = sorted(glob(simset+"pt-*/"))
  tt2  = tqdm(sims,position=1,leave=False)

  mass_arr = []
  sne_arr  = []
  inter_arr = []
  

  for sim_number,sim in enumerate(tt2):
    state_fnames = sorted(glob(sim+"*-state-*.zst"))
    yieldsname = glob(sim + "/*yields.ubj.zst")[0]
    first_state_fname = state_fnames[0]
    last_state_fname = state_fnames[-1]
    first_state = read_state(first_state_fname)
    final_state = read_state(last_state_fname)
    yields = read_yields(yieldsname)
    final_yields = yields.local_26al_final
    
    initial_cluster = first_state.cluster
    final_cluster  = final_state.cluster

    metadata = first_state.metadata

    nmassive = len(initial_cluster.mass[initial_cluster.mass >= 13.0 | units.MSun])
    nsne     = len(final_cluster.kicked[final_cluster.kicked == True])

    # Calculate values for stars which received no enrichment

    hm_i,lm_i = get_high_mass_star_indices(initial_cluster)
    nnointer = 0
    ntotal   = 0
    for star_n,y in enumerate(final_yields):
        if initial_cluster[star_n].mass >= 0.1 | units.MSun:
          if initial_cluster[star_n].mass <= 3.0 | units.MSun:
            if y != 0.0:
              nnointer += 1
            ntotal += 1
    
    # nnointer     = len(final_cluster.mass_26al_local[final_cluster.mass_26al_local != 0.0 | units.MSun])
    # ntotal = len(final_cluster.mass)

    mass_arr.append(nmassive)
    sne_arr.append(nsne)
    inter_arr.append(nnointer/ntotal)

  sm = pd.Series(mass_arr)
  ss = pd.Series(sne_arr)
  si = pd.Series(inter_arr)

  filename = metadata.filename

  rc =  metadata.args.rc
  nstar = metadata.args.n

  tt1.write("{} & {:.1f} & {} & ${:.3f}\pm{:.3f}$ & ${:.3f}\pm{:.3f}$ & ${:.3f}\pm{:.3f}$ \\\\".format(filename,rc,nstar,sm.mean(), sm.std()/(sm.count()**0.5),ss.mean(), ss.std()/(sm.count()**0.5),si.mean(), si.std()/(si.count()**0.5)))
