# README

## What is this archive?

## What are the contents?

### `paper_data`
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
    - An additional compression layer would result in significantly longer extraction time for a marginal size decrease 
  - Note: This archive is too large to store on GitHub, as such, please seek out the Sheffield University data repository instead.

### `python_scripts`



## How to use?