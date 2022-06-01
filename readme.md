### vbr_bayesian_inversion

A 3D bayesian inversion for thermodynamic state using the [VBRc](https://github.com/vbr-calc/vbr)

**paper citation**
Paxman, G. J. G., Lau, H. C. P., Austermann, J., Holtzman, B. K., Havlin, C. 2022. Inference of the timescale-dependent apparent viscosity structure in the upper mantle beneath Greenland. Submitted to AGU Advances.

The archived version of this code for the above paper is located at https://doi.org/10.5281/zenodo.6598795.

Further development may occur on Github at https://github.com/guypaxman/vbr_bayesian_inversion.

### Installation

This code requires MATLAB or GNU-Octave and the [VBRc](https://github.com/vbr-calc/vbr). This code also requires some data files not included in the Github repository (see the following section).

### Getting Data

The main input datasets are available in the Zenodo repository. These include .mat files for the GLAD-M25 shear wave velocity model (`Vs_Model.mat`) and the QL6 1D global attenuation model (`QL6_Model.mat`). Also included is the VBRc parameter sweep lookup table (`sweep_box.mat`).

If you do not download the VBRc lookup table, it can be generated via the `make_sweep.m` script before running the main code for the first time (see below).

The seismic tomography data were originally accessed from the IRIS Earth Model Collaboration data products webpages:

GLAD-M25 (Vs): https://doi.org/10.17611/dp/emc.2020.gladm2500.1.
PREM-QL6 (Q): https://doi.org/10.17611/DP/9991844.

### Usage

The main driving scripts are 
* `make_sweep.m`: if you did not download the VBRc lookup table `sweep_box.mat`, you'll need to generate it. 
* `run_with_*_res.m`: these scripts handle loading the observation data, loading the Vs and Q predictions from the VBRc and running the bayesian inversion. Most adjustable parameters are stored in high-level structures in these files. Note that `run_with_original_res.m` will run a high resolution inversion that requires significant RAM.

### Testing

This repository includes a minimal test suite, see `./tests/readme.md` for more details.

  
