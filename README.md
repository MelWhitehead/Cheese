# Cheese
Volcanic Field cheese model code
written by Prof. Mark Bebbington

Accompanying paper: EGUSPHERE-2025-2010 *NHESS Brief Communication: A magma depletion alternative for vent distribution in volcanic fields*

Mark S. Bebbington, Melody G. Whitehead, & Gabor Kereszturi
Volcanic Risk Solutions, School of Agriculture and Environment, Massey University, Palmerston North, 4472, New Zealand
(under review at time of code release)

Abstract:
The location of a volcanic vent controls an eruption’s hazards, intensities, and impact. Current kernel density estimation methods of future vent locations in volcanic fields assume that locations with more past-vents are more likely to produce future-vents. We examine an alternative hypothesis that an eruption depletes the magma source, causing holes or dips in the spatial density estimate for future vent locations. This is illustrated with the Auckland Volcanic Field, Aotearoa-New Zealand, where both magmatic and phreatomagmatic eruptions have occurred, according to the vent location, with the latter resulting in more explosive eruptions and hence hazard.

All code runs in MATLAB as is, or OCTAVE by moving the functions to the top of each file

## Open cheese_master.m

1. Edit the .txt and .dat files to replace coords_vol2.txt and coast.dat with your own volcanic field data (or run as is with provided AVF data)
2. Select which kernel you would like to run via K = # (see options below)
3. Run

### bw_iG_oneout.m
> K == 1: Isotropic Gaussian kernel with no alteration

### bw_disk_oneout.m
> K == 2: Isotropic Gaussian kernel with disk-based depletion

### bw_weight_oneout.m
> K == 3: Inverse-volume weighted isotropic Gaussian kernel 

### bw_r_oneout.m
> K == 4: Distance-weighted Gaussian kernel
