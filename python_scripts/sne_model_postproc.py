"""
Code to calculate 26Al and 60Fe yields from supernovae
"""

import sys
import os
script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")
import al26_plot
from al26_nbody import Yields,State,Metadata
from amuse.lab import SeBa
import matplotlib.pyplot as plt
from amuse.units import units
from glob import glob
import numpy as np
import pandas as pd
from tqdm import tqdm

# First, read in simulation names
simulation_dirnames = sorted(glob("pt-*/pt-*-*-fr/pt-*-*-fr-*/"))

# Very sorry about this one, create a dictionary of lists for all data being written
# This really sucks, I know, but it works, and is relatively fast for work of indeterminate size

proc_data = {}
proc_data["nstars"]              = []
proc_data["rc"]                  = []
proc_data["sim_number"]          = []
proc_data["star"]                = []
proc_data["initial_mass"]        = []
proc_data["mass"]                = []
proc_data["isotope"]             = []
proc_data["model"]               = []
proc_data["yield_ratio_nodecay"] = []
proc_data["yield_ratio_decay"]   = []

for simulation_dirname in simulation_dirnames:
  state_names = sorted(glob(simulation_dirname+"*state-*.pkl.zst"))
  first_state = al26_plot.read_state(state_names[0])
  last_state  = al26_plot.read_state(state_names[-1])
  first_cluster = first_state.cluster
  last_cluster  = last_state.cluster

  # Calculate which massive stars undergo supernovae (and which ones)
  sn_times,sn_masses,sn_keys = al26_plot.calc_sn_times(first_cluster,return_keys=True)
  print(sn_times,sn_keys)
  
  # Because supernovae keys are sorted in the order in which they occur
  sn_processed = 0

  sne_abs_nodiskdecay_26al_tot = np.zeros(len(first_cluster))
  sne_abs_nodiskdecay_60fe_tot = np.zeros(len(first_cluster))
  sne_abs_diskdecay_26al_tot   = np.zeros(len(first_cluster))
  sne_abs_diskdecay_60fe_tot   = np.zeros(len(first_cluster))

  for i in tqdm(range(0,len(state_names)-1)):
    curr_state_name = state_names[i]
    next_state_name = state_names[i+1]
    # Read in first and next timesteps
    curr_state = al26_plot.read_state(curr_state_name)
    next_state = al26_plot.read_state(next_state_name)
    curr_cluster = curr_state.cluster
    # Quickly calculate the anticipated timestep
    t_curr = curr_state.metadata.time
    t_next = next_state.metadata.time
    dt = t_next - t_curr

    hm_id,lm_id = al26_plot.get_high_mass_star_indices(curr_cluster)

    if sn_processed < len(sn_times):
      if sn_times[sn_processed] >= t_curr.value_in(units.myr):
        if sn_times[sn_processed] <= t_next.value_in(units.myr):
          sn_key = sn_keys[sn_processed]
          # Unfortunately method relies on some key matching stuff, which is not very efficient
          for ii in hm_id:
            if curr_cluster[ii].key == sn_key:
              # Get supernoave yield for star
              al26_sn = curr_cluster[ii].sn_yield_26al
              fe60_sn = curr_cluster[ii].sn_yield_60fe
              for jj in lm_id:
                d        = al26_plot.calc_star_distance(curr_cluster[ii],curr_cluster[jj])
                r_disk   = curr_cluster[jj].r_disk
                eta_disk = al26_plot.calc_eta_disk_sne(r_disk,d)
                al26_inj = al26_sn * eta_disk
                fe60_inj = fe60_sn * eta_disk
                # Store as appropriate
                sne_abs_nodiskdecay_26al_tot[jj] += al26_inj.value_in(units.MSun)
                sne_abs_nodiskdecay_60fe_tot[jj] += fe60_inj.value_in(units.MSun)
          sn_processed += 1

    # Now calculate fraction of total lost to radioactive decay for global model
    half_life_26al = 0.717 | units.Myr
    half_life_60fe = 2.600 | units.Myr
    ln2 = 0.693147
    decay_frac_26al = np.exp((-dt * ln2) / half_life_26al)
    decay_frac_60fe = np.exp((-dt * ln2) / half_life_60fe)
    sne_abs_nodiskdecay_26al_tot *= decay_frac_26al
    sne_abs_nodiskdecay_60fe_tot *= decay_frac_60fe

    
    for star_no,star in enumerate(curr_state.cluster):
      if star.disk_alive == True:
        if star.tau_disk >= t_next:
          sne_abs_diskdecay_26al_tot[star_no] = sne_abs_nodiskdecay_26al_tot[star_no]
          sne_abs_diskdecay_60fe_tot[star_no] = sne_abs_nodiskdecay_60fe_tot[star_no]

  # Now that simulation has been processed, write data for each star to the dataset
  meta = first_state.metadata
  filename = meta.filename
  sim_number  = meta.filename.split("-")[-1]
  nstars   = meta.args.n
  rc       = meta.args.rc
  for i,star in enumerate(first_cluster):
    initial_mass = star.mass.value_in(units.MSun)
    mass = last_cluster[i].mass.value_in(units.MSun)
    # I am especially sorry for this bit!!! I ran out of coffee
    isotopes = ["26al","60fe"]
    for isotope in isotopes:
      proc_data["nstars"].append(nstars)
      proc_data["rc"].append(rc)
      proc_data["sim_number"].append(sim_number)
      proc_data["star"].append(i)
      proc_data["initial_mass"].append(initial_mass)
      proc_data["mass"].append(mass)
      proc_data["isotope"].append(isotope)
      proc_data["model"].append("sne")
      if isotope == "26al":
        stable_yield = star.mass_27al.value_in(units.MSun)
        unstable_yield_nodecay = sne_abs_nodiskdecay_26al_tot[i]
        unstable_yield_decay   = sne_abs_diskdecay_26al_tot[i]
      if isotope == "60fe":
        stable_yield = star.mass_56fe.value_in(units.MSun)
        unstable_yield_nodecay = sne_abs_nodiskdecay_60fe_tot[i]
        unstable_yield_decay   = sne_abs_diskdecay_60fe_tot[i]
      proc_data["yield_ratio_nodecay"].append(unstable_yield_nodecay / stable_yield)
      proc_data["yield_ratio_decay"].append(unstable_yield_decay / stable_yield)
  # At the end of every simulation processed, save data
  print("Finished {}".format(simulation_dirname))
  df = pd.DataFrame.from_dict(proc_data)
  print(df)
  df.to_pickle("sne-model-all-sim-ratios.pkl.zst")