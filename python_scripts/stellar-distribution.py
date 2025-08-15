import sys

script_dir = "/local/jweatson/26al-nbody/"
sys.path.append(script_dir)
sys.path.append(script_dir+"/plotting/")

from al26_plot import calc_cdf,use_tex
from al26_nbody import generate_masses
import matplotlib.pyplot as plt 

x1,y1 = calc_cdf(generate_masses(100000,min_mass=0.01,max_mass=150.))
x2,y2 = calc_cdf(generate_masses(100000,min_mass=0.1,max_mass=60.))

use_tex()
plt.figure(figsize=(5,3))
plt.plot(x1,y1,label="Original Masch. IMF")
plt.plot(x2,y2,label="Modified Masch. IMF")
plt.yscale("log")
plt.xscale("log")
plt.xlim(0.009,100)
plt.ylim(1e-6,1.1)
plt.legend()
plt.savefig("masses.pdf")