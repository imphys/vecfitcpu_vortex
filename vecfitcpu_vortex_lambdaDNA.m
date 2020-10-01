% This m-file is for single emitter localization using a vectorial
% dipole PSF
close all;
clearvars;

% set parameters
params = set_parameters_vortex_lambdaDNA;

% Parallel computing flag: (0) for-loop (1) parfor-loop
% - (1) install parallel computing toolbox
params.flg_parallel = 1;

%% read data and segment data
datain = '\lambdaDNA\lambdaDNA_10frames';
filestr = ['read' datain '.tif'];

% spot segmentation
[allspots,roixy,framelist,params] = get_segmentation(params,filestr,25);

% predict zernike values based on segmented spot positions
xfov = params.pixelsize/1E3*(roixy(:,1)-params.FOV/2);
yfov = params.pixelsize/1E3*(roixy(:,2)-params.FOV/2);
[RAstig3,RAstig5,RComa3,RComa5,RCurv5,RTrefoil,RCurv6] = loadPertubations('read/nat/natCoefficients_mlambda_oil');
params.natPredictions = get_natPredictions(xfov,yfov,RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6)/1E3*params.lambda;

%% MLE fitting routine
[thetainit] = initialvalues(allspots,params);
[thetastore,mu,dmudtheta,merit,numiters] = localization(allspots,thetainit,params);
theta = thetastore(:,:,end);

% compute CRLB using the estimated parameters
[crlb,rcondstore] = get_fisher_crlb(params,mu,dmudtheta);

% compute fitting errors (chi2, LLR, and SSE)
[fiterror] = get_fiterror(mu,allspots,params);

% outlier removal based on statistics of maximum log-likelihood, number of
% iterations, and on found position
[outliers] = get_outliers(theta,merit,numiters,params);

% convert xy to field positions and azim-pola to half-sphere
[theta] = roi2fov(theta,roixy,params);

%% Plots

% plot estimated orientations
figure
norm = [pi/180 pi/180 1];
ystr = {'azim (deg)' 'pola (deg)' 'g2'};
set(gcf,'Position',[180 730    713   165])
for ii = 1:3
    subplot(1,3,ii)
    hold on
    hist(theta(end-3+ii,~outliers)/norm(ii))
    xlabel(ystr{ii})
    ylabel('occurrence')
end

% plot fitting errors; chi^2, LLR, SSE
figure
xstr = {'\chi^2' 'LLR' 'SSE'};
set(gcf,'Position',[181   478   713   165])
for ii = 1:3
    subplot(1,3,ii)
    hist(fiterror(ii,~outliers))
    xlabel(xstr{ii})
    ylabel('occurrence')
end

% plot psfs; data (left) and fits (right)
dipshow(cat(2,allspots(:,:,:,~outliers),mu(:,:,:,~outliers)),'lin');
diptruesize(2000)
colormap hot
