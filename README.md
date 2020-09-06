# SelfCalGroupFinder

Basic Usage:
kdGroupFinder inputfile zmin zmax frac_area [fluxlim] [color] [wcenvalues 1-6] [Bsat_values 1-4] [wchi_values 1-4]> out

- inputfile: the galaxy file in which to find the groups (format described below)
- zmin: minimum redshift. Mainly for volume-limited samples. For flux-limited, use value below all galaxies in file
- zmax: maximum redshift. .Mainly for volume-limited samples. For flux-limited, use value above all galaxies in file.
- fluxlim - 1/0 to delineate if the input file is flux-limted or volume-limited.
- color - 1/0 will the input file have color information (for each galaxy, 1=quiescent, 0=star-forming).
- wcenvalues - these are the 6 free parameters that goven the weights on the total group luminosity.
