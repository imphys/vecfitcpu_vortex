function params = set_parameters_zstack_bead
% This function sets all parameters for vectorial PSF calculations

% flags
params.debugmode = 0;
params.flg_parallel = 0;
params.flg_nat = 0;

% libraries
addpath('funfit');
addpath('funextra');

% parameters: NA, refractive indices of medium, cover slip, immersion fluid,
% nominal value of immersion fluid refractive index matching objective lens
% design, nominal free working distance (in nm), distance image plane from
% cover slip (in nm), wavelength (in nm), emitter position (in nm) with
% z-position from image plane, spot footprint (in nm), axial range min/max
% (in nm), flag for axial range by z-position in medium or by z-stage
% position, sampling in pupil with, sampling in image plane, sampling in
% axial direction.
%

% fitting parameters
params.Nitermax = 50;
params.tollim = 1e-6;
params.varfit = 0;

% camera offset and gain
params.offset = 104.2;
params.gain = 1/0.52;

% PSF/optical parameters
params.NA = 1.47;
params.refmed = 1.518;
params.refcov = 1.523;
params.refimm = 1.518;
params.refimmnom = params.refcov;
params.fwd = 120e3;
params.depth = 90;
params.zrange = [-2500,2500];
params.zspread = [-2500,2500];
params.ztype =  'stage'; % 'medium'
params.lambda = 596;
params.lambdacentral = 596;
params.lambdaspread = [596,596];
params.xemit = 0.0;
params.yemit = 0.0;
params.zemit = 0.0;
params.Npupil = 32;
params.pixelsize = 65;
params.samplingdistance = params.pixelsize;
params.Mx = 31;
params.My = params.Mx;
params.Mz = 1;

params.xrange = params.pixelsize*params.Mx/2;
params.yrange = params.pixelsize*params.My/2;

% SAF check
if and(params.NA>params.refmed, params.depth<4*params.lambda)
    [zvals,~] = set_saffocus(params);
else
    [zvals,~] = get_rimismatchpars(params);
end
params.zvals = zvals;

% sanity check on position emitter w.r.t. cover slip
if strcmp(params.ztype,'stage')
    zcheck = params.depth+params.zemit;
end
if strcmp(params.ztype,'medium')
    zmin = params.zrange(1);
    zcheck = zmin+params.depth+params.zemit;
end
if (zcheck<0)
    fprintf('Warning! Emitter must be above the cover slip:\n')
    fprintf('Adjust parameter settings for physical results.\n')
end

% sanity check on refractive index values
if (params.NA>params.refimm)
    fprintf('Warning! Refractive index immersion medium cannot be smaller than NA.\n')
end
if (params.NA>params.refcov)
    fprintf('Warning! Refractive index cover slip cannot be smaller than NA.\n')
end

% parameters needed for fixed dipole PSF only: emitter/absorber dipole
% orientation (characterized by angles pola and azim)
params.dipoletype = 'diffusion';
params.pola = 90.0*pi/180;
params.azim = 0.0*pi/180;

% diffusion coefficient
welldepth = eps;
g2 = (3+welldepth^2-3*welldepth*coth(welldepth))/welldepth^2;
params.welldepth = welldepth;
params.g2 = g2;

% aberrations (Zernike orders [n1,m1,A1,n2,m2,A2,...] with n1,n2,... the
% radial orders, m1,m2,... the azimuthal orders, and A1,A2,... the Zernike
% coefficients in lambda rms, so 0.072 means diffraction limit)
% params.aberrationcorrected = false;
% params.aberrationsoffset = [];
params.aberrations = [ ...
    2, -2,  0;
    2,  2,  0;
    3, -1,  0;
    3,  1,  0;
    3, -3,  0;
    3,  3,  0;
    4,  0,  0;
    4, -2,  0;
    4,  2,  0;
    5, -1,  0;
    5,  1,  0;
    6,  0,  0];
params.aberrations(:,3) =  params.aberrations(:,3)*params.lambdacentral;

% DOE/SLM
params.doetype = 'none';
% params.doetype = 'vortex';
% params.ringradius = 1;
% params.doelevels = 64;
% params.zonefunction = params.aberrations;
% params.doephasedepth = 593;
% params.doevortexflip = 1;

% Bead parameters for convolution with PSF and derivatives, beaddiameter in nm
params.bead = true;
params.beaddiameter = 180;
% check on meaningfullness of bead convolution
if params.beaddiameter<params.pixelsize
    params.bead = false;
end

% Fit model parameters: signal photon count, background photons/pixel, read
% noise variance for sCMOS camera's, fit model, output labels depending on
% fit modelnn
params.readnoisestd = 0;
params.readnoisevariance = 0;

% model parameters
params.alpha = 0;
params.beta = 0;
params.K = 1;
params.m = 1;
params.excitation = 'zstack';

% fitting model
params.fitmodel = 'xyz-aberrations';
if strcmp(params.fitmodel,'xy')
    params.numparams = 4;
elseif strcmp(params.fitmodel,'xyz')
    params.numparams = 5;
elseif strcmp(params.fitmodel,'xyz-aberrations')
    params.numparams = 5+numel(params.aberrations(:,3));
elseif strcmp(params.fitmodel,'xy-azim')
    params.numparams = 5;
elseif strcmp(params.fitmodel,'xy-azim-pola')
    params.numparams = 6;
elseif strcmp(params.fitmodel,'xyz-azim-pola')
    params.numparams = 7;
elseif strcmp(params.fitmodel,'xy-azim-diffusion')
    params.numparams = 6;
elseif strcmp(params.fitmodel,'xy-azim-pola-diffusion')
    params.numparams = 7;
elseif strcmp(params.fitmodel,'xyz-azim-pola-diffusion')
    params.numparams = 8;
elseif strcmp(params.fitmodel,'xyz-azim-pola-aberrations')
    params.numparams = 7+numel(params.aberrations(:,3));
elseif strcmp(params.fitmodel,'xyz-azim-pola-diffusion-aberrations')
    params.numparams = 8+numel(params.aberrations(:,3));
end


% % calculate auxiliary vectors for chirpz
PupilSize = params.NA/params.lambda;
[Ax,Bx,Dx] = prechirpz(PupilSize,params.xrange,params.Npupil,params.Mx);
[Ay,By,Dy] = prechirpz(PupilSize,params.yrange,params.Npupil,params.My);

params.Axmt = repmat(Ax,params.Mx,1);
params.Bxmt = repmat(Bx,params.Mx,1);
params.Dxmt = repmat(Dx,params.Mx,1);
params.Aymt = repmat(Ay,params.Npupil,1);
params.Bymt = repmat(By,params.Npupil,1);
params.Dymt = repmat(Dy,params.Npupil,1);

params.cztN = params.Npupil;
params.cztM = params.Mx;
params.cztL = params.Npupil+params.Mx-1;

% show intermediate results for monitoring code
params.debugmode = 0;