function parameters = set_parameters_xy
% This function sets all parameters for vectorial PSF calculations
%
% copyright Sjoerd Stallinga, TU Delft, 2017
%
% parameters: NA, refractive indices of medium, cover slip, immersion fluid,
% wavelength (in nm), nominal emitter position (in nm) with z-position from
% cover slip-medium interface, spot footprint (in nm), axial range (in nm),
% sampling in pupil with (must be even), sampling in image plane (must be odd),
% sampling in axial direction, emission or excitation light path.

parameters.NA = 1.33; % Numerical Aperture objective
parameters.refmed = 1.33; % refractive index medium
parameters.refcov = 1.52; % refractive index cover slip
parameters.refimm = 1.52; % refractive index immersion fluid
parameters.refimmnom = parameters.refcov; % nominal refraction index immersion fluid
parameters.fwd = 140e3; % free working distance objective
parameters.depth = 0; % depth inside medium for nominal focus 
parameters.zstage = 0.0; % z-stage position
parameters.zrange = [0,0]; % range of z-values in z-estimation
parameters.zspread = [0,0]; % range of ground truth z-values in z-estimation
parameters.ztype = 'medium'; % z-position estimated in medium
% parameters.ztype = 'stage'; % z-position estimated as z-stage position
parameters.lambda = 550.0; % emission wavelength (nm)
parameters.lambdacentral = 550; % nominal wavelength for DOEs
parameters.lambdaspread = [550,550]; % range of ground truth wavelengths in wavelength estimation
parameters.xemit = 0.0; % x-position emitter
parameters.yemit = 0.0; % y-position emitter
parameters.zemit = 0.0; % z-position emitter
parameters.Npupil = 56; % # pupil sampling points, now L=N+M-1=64 which enables power2 fft
parameters.Npupil = 40; % # pupil sampling points, now L=N+M-1=48=3*16 which enables power2,3 fft
% parameters.Npupil = 24; % now L=N+M-1=32 which enables power2 fft
parameters.fitmode = 'OTF'; % use fast OTF based fitting, only if parameters.fitmodel = 'xy' !!
% parameters.fitmode = 'PSF'; % or direct vector PSF model
switch parameters.fitmode
  case 'OTF'
    parameters.psfsamplingdistance = 0.25*(parameters.lambda/parameters.NA/4); % sampling distance in PSF space
    parameters.Mpsfx = 129; % #sampling points in PSF space
    parameters.Mpsfy = 129; % #sampling points in PSF space
    parameters.psfxrange = parameters.psfsamplingdistance*parameters.Mpsfx/2; % 1/2 of size PSF space in x
    parameters.psfyrange = parameters.psfsamplingdistance*parameters.Mpsfy/2; % 1/2 of size PSF space in y
    parameters.Notf = 40; % #sampling points in OTF space
    parameters.otfsamplingdistance = 4*parameters.NA/parameters.lambda/parameters.Notf; % sampling distance in OTF space
    parameters.otfxrange = parameters.otfsamplingdistance*parameters.Notf/2; % 1/2 of size OTF space in x
    parameters.otfyrange = parameters.otfsamplingdistance*parameters.Notf/2; % 1/2 of size OTF space in y
    parameters.pixelsize = 80; % camera pixel size
    parameters.roisamplingdistance = parameters.pixelsize; % sampling distance in ROI
    parameters.Mroix = 9; % #pixels in ROI in x
    parameters.Mroiy = 9; % #pixels in ROI in y
    parameters.xroirange = parameters.roisamplingdistance*parameters.Mroix/2; % 1/2 of size PSF space in x
    parameters.yroirange = parameters.roisamplingdistance*parameters.Mroiy/2; % 1/2 of size PSF space in y
    
    parameters.samplingdistance = parameters.roisamplingdistance; % sampling distance in ROI
    parameters.Mx = parameters.Mroix;
    parameters.My = parameters.Mroiy;
  case 'PSF'
    parameters.pixelsize = 80; % camera pixel size
    parameters.samplingdistance = parameters.pixelsize; % sampling distance in ROI
    parameters.Mx = 9; % #pixels in ROI in x
    parameters.My = 9; % #pixels in ROI in y
end
parameters.xrange = parameters.samplingdistance*parameters.Mx/2; % 1/2 of ROI size in x
parameters.yrange = parameters.samplingdistance*parameters.My/2; % 1/2 of ROI size in y
parameters.Mz = 1; % #focal slices (=1 in xy-mode)
parameters.Mlambda = 1; % #wavelengths (=1 in xy-mode)

% sanity check on position emitter w.r.t. cover slip
if strcmp(parameters.ztype,'stage')
  zcheck = parameters.depth+parameters.zemit;
end
if strcmp(parameters.ztype,'medium')
  zmin = parameters.zrange(1);
  zcheck = zmin+parameters.depth+parameters.zemit;
end
if (zcheck<0)
  fprintf('Warning! Emitter must be above the cover slip:\n')
  fprintf('Adjust parameter settings for physical results.\n')  
end

% sanity check on refractive index values
if (parameters.NA>parameters.refimm)
  fprintf('Warning! Refractive index immersion medium cannot be smaller than NA.\n')
end
if (parameters.NA>parameters.refcov)
  fprintf('Warning! Refractive index cover slip cannot be smaller than NA.\n')
end

% parameters needed for fixed dipole PSF only: emitter/absorber dipole
% orientation (characterized by angles pola and azim), detection/illumination
% polarization in objective lens back aperture (characterized by angles
% alpha and beta).
parameters.dipoletype = 'free';
% parameters.dipoletype = 'fixed';
parameters.pola = 45.0*pi/180;
parameters.azim = 0.0*pi/180;
parameters.polarizationpupil = false;
parameters.alpha = 45.0*pi/180;
parameters.beta = 45.0*pi/180;

% aberrations (Zernike orders [n1,m1,A1,n2,m2,A2,...] with n1,n2,... the
% radial orders, m1,m2,... the azimuthal orders, and A1,A2,... the Zernike
% coefficients in lambda rms, so 0.072 means diffraction limit)
parameters.aberrations = [2,0,0.0; 2,2,0.0; 2,-2,0.0; 3,1,0.0; 3,-1,0.0; 4,0,0.0];
parameters.aberrations(:,3) =  parameters.aberrations(:,3)*parameters.lambdacentral;
parameters.aberrationoffset = false;
parameters.aberrations_delta = parameters.aberrations; 
parameters.aberrations_delta(:,3) = -parameters.aberrations(:,3);

% Diffractive Optical Element (DOE) parameters: type (binary, blaze,
% sinusiodal, interpolant), phase depth, zone function (as Zernike list)

parameters.doetype = 'none';
% parameters.doetype = 'binary';
% parameters.doetype = 'blazed';
% parameters.doetype = 'sinusoidal';
% parameters.doetype = 'spindle';
% parameters.doetype = 'azimuthramp';
% parameters.doephasedepth = 0.5*parameters.lambdacentral;
% parameters.interpolant = [0.0 0.1 0.2 0.3]*parameters.lambdacentral;
% numsingularpoints = 9;
% singularpointspacing = 0.633;
% singlength = singularpointspacing*(numsingularpoints-1)/2;
% parameters.singularpoints = linspace(-singlength,singlength,numsingularpoints);
% numzones = 2;
% parameters.ringradius = (linspace(0,1,numzones+1)).^(1/2);
% parameters.ringradius = [0.0,0.4627,0.9044,1.0];
% parameters.ringradius = [0.0,0.8887,1.0];
% parameters.zonefunction = [0,0,0.0; 1,1,0.0; 1,-1,0.7; 2,0,0.0; 4,0,-0.0; 2,-2,0; 4,-2,-0.0];

% Fit model parameters: signal photon count, background photons/pixel, read
% noise variance for sCMOS camera's, fit model, output labels depending on
% fit model

parameters.signalphotoncount = 2000;
parameters.backgroundphotoncount = 5;
parameters.readnoisestd = 0;

parameters.fitmodel = 'xy';
parameters.numparams = 4;
parameters.outputlabels = {'x','y','N_p_h','b_g'};
parameters.outputunits = {' nm',' nm',' ',' '};
parameters.weightvector = [1 1 0 0];
parameters.weightvector = parameters.weightvector/sum(parameters.weightvector);

% Bead parameters for convolution with PSF and derivatives, beaddiameter in nm
parameters.bead = false;
parameters.beaddiameter = 200;

% check on meaningfullness of bead convolution
if parameters.beaddiameter<=parameters.pixelsize
  parameters.bead = false;
end

% show intermediate results for monitoring code
parameters.debugmode = 0;

end

