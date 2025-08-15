import sys
import os
script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")
import al26_plot
from al26_nbody import Yields,State,Metadata
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import pandas as pd
from tqdm import tqdm
from amuse.units import units

dirii = glob("pt-*/pt-*-*-fr/pt-*-*-fr-*/")
rc_arr = ["0.3","1.0","3.0"]
ns_arr = ["100","1000","10000"]
markers = ["o","v","s"]

al26_plot.use_tex()
plt.figure(figsize=(7.5,3))

for i,rc in enumerate(rc_arr):
  for j,ns in enumerate(ns_arr):
    dirii = glob("pt-{}/pt-{}-{}-fr/pt-*/".format(rc,rc,ns))
    zen_arr = []
    mrho_arr = []
    for k,diri in enumerate(dirii):
      print(i,j,k)
      statename  = glob(diri + "/pt-*-*-fr-*-state-00000.pkl.zst")[0]
      yieldsname = glob(diri + "/*yields.ubj.zst")[0]
      state = al26_plot.read_state(statename)
      cluster = state.cluster 
      local_densities = al26_plot.calc_local_densities(cluster)
      med_local_dens  = np.median(local_densities)
      yields = al26_plot.read_yields(yieldsname)
      final_yields = yields.local_26al_final

      # print(yields)

      # print(cluster.mass.value_in(units.MSun))
      # masses = cluster.mass.value_in(units.MSun)
      # print(np.argwhere((masses >= 0.1) | (masses <= 3.0)))

      nstars = 0
      ninter = 0
      for star_n,y in enumerate(final_yields):
        if cluster[star_n].mass >= 0.1 | units.MSun:
          if cluster[star_n].mass <= 3.0 | units.MSun:
            if y != 0.0:
              ninter += 1
            nstars += 1
      zen = ninter / nstars
      zen_arr.append(zen)
      mrho_arr.append(med_local_dens)
    label  = "$r_c$="+rc+", $N_\\star$="+ns
    plt.scatter(mrho_arr,zen_arr,color="C"+str(i),marker=markers[j],label=label,s=3)

plt.xscale("log")
plt.ylabel("Disk enrichment fraction, $Z_\\mathrm{local,enrich}$")
plt.xlabel("Region median local density, $\\tilde{\\rho}$")
plt.grid(ls=":")
plt.ylim(0,1)
plt.legend(loc="upper center",ncols=3,prop={'size': 6})
plt.savefig("median-local-density-versus-disk-enrich-local.pdf",bbox_inches="tight")