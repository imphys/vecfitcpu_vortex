# Software for Simultaneous orientation and 3D localization microscopy with a Vortex point spread function

This code is distributed as accompanying software for the article *Simultaneous orientation and 3D localization microscopy with a Vortex point spread function* by Christiaan N. Hulleman, Rasmus Ø. Thorsen, Sjoerd Stallinga, and Bernd Rieger.

## General concept

The code uses a vectorial point-spread-function model to perform a maximum likelihood estimate on single-molecule emitters. The parameters of interest are the emitter position, photon counts, and orientation together with rotational constraint. The found values are used to calculate the Cramer-Rao Lower Bound (CRLB) for each parameter, and the CRLBs are returned along with the estimated parameters.

## Example Usage
Several examples of how to use the code on simulated and experimental data are shown in the following MATLAB scripts

- vecfitcpu_vortex_simfits.m
- vecfitcpu_vortex_lambdaDNA.m
- vecfitcpu_zstack_bead.m
- generate_zernike_surfaces.m

Additional experimental data can be found via https://doi.org/10.4121/c.5136125.

## MATLAB
The code is written in MATLAB, and tested to work in MATLAB R2018-R2020. The DIPImage toolbox for MATLAB is required, please see http://www.diplib.org for installation instructions.


### Further questions
For further questions feel free to create an issue on GitHub. You can also contact the authors:
