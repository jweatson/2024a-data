import sys
import os

script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")

import matplotlib.colors as mcolors
import matplotlib
from al26_plot import read_yields,read_state,calc_cdf,use_tex,calc_sn_times
from al26_nbody import Yields,State,Metadata,generate_masses
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
from amuse.lab import SeBa
from amuse.datamodel import Particles


simsets = sorted(glob("./pt-*/pt*/"))

masses = []
taus   = []

for simset in simsets:
  sims = sorted(glob(simset+"pt-*/"))
  max_al_local = []
  
  for sim_number,sim in enumerate(tqdm(sims)):
    yields_fname = sorted(glob(sim+"*yields*.zst"))[-1]
    state_fnames = sorted(glob(sim+"*-state-*.zst"))
    first_state_fname = state_fnames[0]
    first_state = read_state(first_state_fname)
    metadata = first_state.metadata
    nstars = metadata.args.n
    cluster = first_state.cluster

    hm_cluster = cluster[cluster.mass > 13.0 | units.MSun]
    lm_cluster = cluster[cluster.mass <= 3.0 | units.MSun]
    for star in hm_cluster:
      masses.append(star.mass.value_in(units.MSun))
    for star in lm_cluster:
      taus.append(star.tau_disk.value_in(units.Myr))

masses = np.asarray(masses)
taus   = np.asarray(taus)

stellar = SeBa(number_of_workers=8)
# detect_supernova = stellar.stopping_conditions.supernova_detection
# detect_supernova.enable()

nstars = len(masses)
ndisks = len(taus)
stars = Particles(nstars)
stars.mass = masses | units.MSun

print("Adding {} particles".format(nstars))
stellar.particles.add_particle(stars)
print("Done!")

t = 0.0 | units.myr

print(stellar.particles)

ages = np.linspace(0,20,256) | units.myr

frac_sn = []
frac_sn_viable = []

from amuse.datamodel import Particles
limit = Particles(1)
limit.mass = [25.0] | units.MSun
sne_yield_limit = calc_sn_times(limit)
print(sne_yield_limit)


for age in tqdm(ages):
  stellar.evolve_model(age)
  sne = len(stellar.particles[stellar.particles.stellar_type >= 13 | units.stellar_type])
  frac_sn.append(1-(sne / nstars))

frac_tau = []
for age in tqdm(ages):
  tfrac = (taus > age.value_in(units.myr)).sum() / ndisks
  frac_tau.append(tfrac)

def quick_interp(x,y,xx):
  from scipy.interpolate import Akima1DInterpolator
  interp = Akima1DInterpolator(x,y)
  yy = float(interp(xx))
  return yy 

hmm = quick_interp(ages.value_in(units.myr),np.asarray(frac_sn),sne_yield_limit[0])
dmm = quick_interp(ages.value_in(units.myr),np.asarray(frac_tau),sne_yield_limit[0])

print("HM fraction remaining at 25msol limit: {}".format(hmm))
print("Disk fraction remaining after 25msol limit: {}".format(dmm))
print("Done, plotting")

use_tex()
fig,ax1 = plt.subplots(figsize=(7.5,3))
ax2 = plt.twinx()

ax1.plot(ages.value_in(units.myr),frac_tau,"tab:blue")
ax1.tick_params(axis="y",colors="tab:blue")
ax1.set_ylabel("Remaining disk fraction, $Z_\\mathrm{disk}$",color="tab:blue")

ax1.axvline(sne_yield_limit[0],color="k",ls=":")
ax1.axhline(dmm,color="tab:blue",ls=":")
ax2.axhline(hmm,color="tab:orange",ls=":")

ax2.plot(ages.value_in(units.myr),frac_sn,"tab:orange")
ax2.tick_params(axis="y",colors="tab:orange")
ax2.set_ylabel("Reamining HM star fraction, $Z_{\\star,\\mathrm{HM}}$",color="tab:orange")

ax1.set_xlabel("Time, $t$ (Myr)")

ax1.set_xlim(0,20)
ax1.set_ylim(-0.05,1.05)
ax2.set_ylim(-0.05,1.05)
ax1.grid(ls=":")

plt.savefig("SNe-fraction.pdf",bbox_inches="tight")
