# SelfCalGroupFinder

Basic Usage:
kdGroupFinder inputfile zmin zmax frac_area [fluxlim] [color] [wcenvalues 1-6] [Bsat_values 1-4] [wchi_values 1-4]> out

- inputfile: the galaxy file in which to find the groups (format described below)
- zmin: minimum redshift. Mainly for volume-limited samples. For flux-limited, use value below all galaxies in file
- zmax: maximum redshift. .Mainly for volume-limited samples. For flux-limited, use value above all galaxies in file.
- fluxlim - 1/0 to delineate if the input file is flux-limted or volume-limited.
- color - 1/0 will the input file have color information (for each galaxy, 1=quiescent, 0=star-forming).
- wcenvalues - these are the 6 free parameters that govern the weights on the total group luminosity. These are taken from Equation (4) in the SDSS Application paper (Tinker 2020).
 - wcen[1] - omega_L,sf 
 - wcen[2] - sigma_sf
 - wcen[3] - omega_L,q
 - wcen[4] - sigma_q
 - wcen[5] - omega_0,sf
 - wcen[6] - omega_0,q
- Bsat_values - these are the 4 free parameters that set the satellite probability threshold. Taken from Equation (3).
 - Bsat[1] - beta_0,q
 - Bsat[2] - beta_L,q
 - Bsat[3] - beta_0,sf
 - Bsat[4] - beta_L,sf
- wchi values - these are the 4 free parameters that set the weights by normalized galaxy property chi. Taken from Equation (5).
  - wchi[1] - omega_chi,0,sf
  - wchi[2] - omega_chi,0,q
  - wchi[3] - omega_chi,L,sf
  - wchi[4] - omega_chi,L,q
