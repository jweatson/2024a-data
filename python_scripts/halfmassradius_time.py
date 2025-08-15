import sys
import os
script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")
import al26_plot
from al26_nbody import Yields,State,Metadata,myr
from amuse.units import units
import matplotlib.pyplot as plt
from matplotlib import cm
from glob import glob
import numpy as np
import pandas as pd
from tqdm import tqdm

rc_arr = ["1.0"]
ns_arr = ["100"]

al26_plot.use_tex()
plt.figure(figsize=(7.5,3))

dirii = glob("pt-*/pt-*-*-fr/pt-*-*-fr-*/")
for i,rc in enumerate(rc_arr):
  for j,ns in enumerate(ns_arr):
    dirii = glob("pt-{}/pt-{}-{}-fr/pt-*/".format(rc,rc,ns))
    for k,diri in enumerate(tqdm(dirii)):
      statenames = sorted(glob(diri + "/pt-*-*-fr-*-state-*.pkl.zst"))
      t_arr   = []
      hmr_arr = []
      cmap    = cm.winter(np.linspace(0,1,len(dirii)))
      for l,statename in enumerate(statenames):
        state = al26_plot.read_state(statename)
        time    = state.metadata.time.value_in(myr)
        cluster = state.cluster
        # hmr = al26_plot.calc_cluster_half_mass(cluster)
        hmr = cluster.virial_radius().value_in(units.pc)
        t_arr.append(time)
        hmr_arr.append(hmr)
      plt.plot(t_arr,hmr_arr,color=cmap[k])
        
plt.xlim(0,20)
plt.ylim(0,3.5)
plt.grid(ls=":")
plt.xlabel("Simulation time, $t$ (Myr)")
plt.ylabel("Region half-mass radius, $r_{1/2}$ (pc)")

plt.savefig("hmr-vs-time-100-1.0-vir.pdf",bbox_inches="tight")