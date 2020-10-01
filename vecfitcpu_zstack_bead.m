% This m-file is for single emitter localization using a vectorial
% dipole PSF
close all;
clearvars;

% set parameters
params = set_parameters_zstack_bead;

%% read data (segmented beda)
file = '\bead\bead180_oil';
datain = ['read' file '.tif'];

InfoImage =  imfinfo(datain);
try
    allspots = double(LoadTiff16bit(datain, [1 length(InfoImage)]));
catch
    allspots = zeros(InfoImage(1).Width,InfoImage(1).Height,length(InfoImage));
    for j=1:length(InfoImage)
        allspots (:,:,j) = double(imread(datain, j));
    end
    disp('Warning: Used slower imread instead of LoadTiff')
end
allspots = (allspots-params.offset)/params.gain;
allspots(allspots<=0) = 1e-3;

% update parameters
params.K = size(allspots,3); % number of z-slices
params.Ncfg = size(allspots,4); % number of configurations

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

%% Plots

% plot zernike values
figure
set(gcf,'Position',[120 700 555 250])
hold on; box on;
numzers = params.numparams-5;
plot(0:numzers+1,zeros(1,numzers+2),'-','Color',[.85 .85 .85],'LineWidth',0.5)
orders = params.aberrations(:,1:2);
allxticks = 1:numzers;
allxticklabels = cell(numzers,1);
for jzer = 1:numzers
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end
plot(1:numzers,theta(6:end,:)/params.lambda*1E3,'k-*','MarkerSize',5)
xticks(allxticks)
xtickangle(25)
xticklabels(allxticklabels)
xlim([0 numzers+1])
xlabel('zernike mode (n,m)');
ylabel('rms value (m\lambda)');

% plot psfs; bead (left) and fits (right)
dipshow(cat(2,allspots,mu),'lin');
diptruesize(1000)
colormap hot
