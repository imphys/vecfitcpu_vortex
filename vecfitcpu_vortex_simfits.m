% This m-file is for single emitter localization using a vectorial
% dipole PSF
close all;
clearvars;

% set parameters
params = set_parameters_vortex_sim;
params.Ncfg = 500; % generate and fit Ncfg randomized model PSFs.

% Parallel computing flag: (0) for-loop (1) parfor-loop
% - (1) install parallel computing toolbox otherwise use (0)
params.flg_parallel = 1;

%% generate ground truth PSFs
% set pupil function
[wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);

allpsfs = zeros(params.Mx,params.My,params.K,params.Ncfg);
object = zeros(params.numparams,params.Ncfg);
for jj = 1:params.Ncfg
    
    % true parameters
    dx = (1-2*rand)*params.pixelsize;
    dy = (1-2*rand)*params.pixelsize;
    dz = 0;
    Nphotons = 4000;
    Nbackground = 10;
    dazim = pi*rand;
    dpola = acos(1-2*rand);
    dg2 = 0.75;
    
    % generate object and PSFs
    object(:,jj) = [dx dy dz Nphotons Nbackground dazim dpola dg2];
    allpsfs(:,:,params.K,jj) = poissonrate(params,object(:,jj),PupilMatrix,allzernikes,wavevector,wavevectorzimm);
    
end

% add noise
allspots = 1e12*imnoise(allpsfs*1e-12,'poisson');

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

% compute the mean and standard deviation of the results
[thetafinal,thetamean,thetastd,crlbmean] = get_statistics(params,object,theta,crlb,outliers);

%% Plot results

% plot parameter errors; bias +- std (red) and CRLB (black)
figure
norm = [1 1 1 1 1 pi/180 pi/180 1];
ystr = {'dx (nm)' 'dy (nm)' 'dz (nm)' 'dN (photons)' 'db (photons/pixel)' 'dazim (deg)' 'dpola (deg)' 'dg2'};
set(gcf,'Position',[180 730 1615 250])
for ii = 1:params.numparams
    subplot(1,params.numparams,ii)
    hold on
    plot(thetafinal(ii,:)/norm(ii),'.','Color',[1 1 1]*0.75,'LineWidth',1.5)
    pmean=plot(thetamean(ii)*ones(params.Ncfg,1)/norm(ii),'r','LineWidth',1.5);
    plot((thetamean(ii)+thetastd(ii))*ones(params.Ncfg,1)/norm(ii),'r','LineWidth',1.5)
    plot((thetamean(ii)-thetastd(ii))*ones(params.Ncfg,1)/norm(ii),'r','LineWidth',1.5)
    pcrlb=plot((thetamean(ii)+crlbmean(ii))*ones(params.Ncfg,1)/norm(ii),'-.k','LineWidth',1.5);
    plot((thetamean(ii)-crlbmean(ii))*ones(params.Ncfg,1)/norm(ii),'-.k','LineWidth',1.5)
    xlabel('cfg')
    ylabel(ystr{ii})
end
legend([pmean,pcrlb],'Mean +- std','Mean +- CRLB')

% plot fitting errors
figure
xstr = {'\chi^2' 'LLR' 'SSE'};
set(gcf,'Position',[181   478   713   165])
for ii = 1:3
    subplot(1,3,ii)
    hist(fiterror(ii,:))
    xlabel(xstr{ii})
    ylabel('occurrence')
end

% plot psfs; ground truth (left) and fits (right)
dipshow(cat(2,allspots,mu),'lin');
diptruesize(2000)
colormap hot

