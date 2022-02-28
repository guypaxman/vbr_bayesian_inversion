## vbr_bayesian_inversion

A 3D bayesian inversion for thermodynamic state using the [VBRc](https://github.com/vbr-calc/vbr)


**paper citation and link here**

The archived version of this code for the above papers is at **ZENODO DOI**. 

Further development may occur on github at https://github.com/guypaxman/vbr_bayesian_inversion

## Installation

This code requires MATLAB or GNU-Octave and the [VBRc](https://github.com/vbr-calc/vbr). This code also requires some datafiles not included in the github repository (see the following section).

### Getting Data 

Info here on downloading the .mat GLAD25 data. The VBRc lookup table can also be downloaded or generated via the `make_sweep.m` script.

### Usage

The main driving scripts are 
* `make_sweep.m` : if you did not download the VBRC lookup table `sweep_box.mat`, you'll need to generate it. 
* `run_with_*_res.m` : these scripts handle loading the observation data, loading the Vs and Q predictions from the VBRc and running the bayesian inversion. Most adjustable parameters are stored in high-level structures in these files. Note that `run_with_original_res.m` will run a high resolution inversion that requires significant RAM.

### Testing

This repository includes a minimal test suite, see `./tests/readme.md` for more details.

  
