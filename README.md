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

Example Usage, taken randomly from the posterior distribution of parameters in the MCMC anlaysis in the paper.
- kdGroupFinder_omp sdss_fluxlim_v1.0.dat  0 1 0.179 1 1  13.1 2.421701964499512 12.9 4.84 17.4 2.67 -0.92 10.25 12.993 -8.04 2.68 1.10 2.23 0.48 >outxx

NB: The code expects a tabulated halo mass function is in the run directory, in a file called "halo_mass_function.dat." I have supplied one in the repo for the Bolshoi Planck cosmology using the Tinker et al (2008) halo mass function.

NB: The code will run for 5 iterations and then output the current state of the group catalog. Five is usually a reasonable number for convergence of the satellite fraction to a couple of percent. This can be modified inside the code. There are several features currentl "turned off," which involve populating the halos of a simulation with HODs, and tabulating the L_sat values. These can be turned on with little effort in the main() function, but user-defined input is required to make them actually run, by supplying the necessary files.
