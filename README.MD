# README

This is an archive of the data and plotting scripts used in the paper _Towards a unified injection model of short-lived radioisotopes in N-body simulations of star-forming regions_, I am coming up to the end of this postdoc, and just wanted to make sure it is somewhere.

A more complete data set has been backed up privately, and is available upon request. 

# `python_scripts`

A collection of scripts that were used to plot the data in the paper. Some of these are quite sloppily written, so your mileage as to whether it works on your machine varies.

# `paper_data`
- `data.tar`'
  - Contains all simulation output data:
    - Simulations are indexed in the form:
      - `pt-<region_radius>/pt-<region_radius>-<star_count>-fr/pt-<region_radius>-<star_count>-fr-<sim_number>`
    - This tarball contains the following common files
      - `state` files: contains all star positions and disk properties, a snapshot of the simulation at a given time
        - Also contains simulation metadata
        - 100 saved per simulation
        - Stored as a ZSTD compressed Pickle file, can be read in using the functions in al26_plot
      - `yield` files: contains final and time-based SLR yields for each star in the simulation, one per simulation
      - `cluster-yield.csv`: total yields over time, used to show that simulation is operating as it is running
  - This has been tarballed due to the extremely large number of individual files (~30,000 files)
    - No compression has been applied, as all of the data has already been compressed with ZSTD
    - An additional compression layer would result in significantly longer extraction time for a marginal size decrease.



## How to use?