# Software for Simultaneous orientation and 3D localization microscopy with a Vortex point spread function

This code is distributed as accompanying software for the article *Simultaneous orientation and 3D localization microscopy with a Vortex point spread function* by Christiaan N. Hulleman, Rasmus Ø. Thorsen, Sjoerd Stallinga, and Bernd Rieger.

Any reuse of this code should cite the original associated publication. 

## General concept

The code uses a vectorial point-spread-function model to perform a maximum likelihood estimate on single-molecule emitters. The parameters of interest are the emitter position, photon counts, and orientation together with rotational constraint. The found values are used to calculate the Cramer-Rao Lower Bound (CRLB) for each parameter, and the CRLBs are returned along with the estimated parameters.

## Using vecfitcpu
Several examples of how to use the code on simulated and experimental data are shown in the following MATLAB scripts

- vecfitcpu_vortex_simfits.m
- vecfitcpu_vortex_lambdaDNA.m
- vecfitcpu_zstack_bead.m
- generate_zernike_surfaces.m

Additional experimental data can be found via https://doi.org/10.4121/c.5136125.

## MATLAB
The code is written in MATLAB, and tested to work in MATLAB R2018-R2020. The DIPImage toolbox for MATLAB is required, please see http://www.diplib.org for installation instructions.


## Further questions
For further questions feel free to create an issue on GitHub. You can also contact the authors:

Christiaan N. Hulleman (c.n.hulleman@tudelft.nl)

Rasmus Ø. Thorsen (r.o.thorsen@tudelft.nl)

Sjoerd Stallinga (s.stallinga@tudelft.nl)

Bernd Rieger (b.rieger@tudelft.nl)

## Reference

If you find this code useful for your research, please cite
```
@article {Hulleman2020.10.01.322834,
	author = {Hulleman, Christiaan N. and Thorsen, Rasmus {\O}. and Stallinga, Sjoerd and Rieger, Bernd},
	title = {Simultaneous orientation and 3D localization microscopy with a Vortex point spread function},
	year = {2020},
	doi = {10.1101/2020.10.01.322834},
	journal = {bioRxiv}
}
```
