% This m-file is for field-dependent aberrations
close all;clearvars;

% load bead fits
load('read/bead/bead180_fits')
addpath('funnat');

fov = params.FOV;
pixelsize = params.pixelsize;

% estimated zernike amplitudes
zers = theta(6:end,:)';

% generate grid
xn = pixelsize/1E3*theta(1,:)'; % (um)
yn = pixelsize/1E3*theta(2,:)'; % (um)
xim = pixelsize/1E3*((1:fov)-fov/2);
yim = pixelsize/1E3*((1:fov)-fov/2);
xim = xim(1:10:end); % downsample grid
yim = yim(1:10:end);
[Xim,Yim] = meshgrid(xim,yim);

% get nat coefficients, surfaces and predictions
[RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6] = get_natCoefficients(xn,yn,zers,params);
[Zsurface] = get_natSurfaces(Xim,Yim,RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6);
[Zpredict] = get_natPredictions(xn,yn,RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6);

%% Plot

% plot estimated zernike coefficients
figure
set(gcf,'Position',[1285 720 555 250])
hold on; box on;
numzers = params.numparams-5;
plot(0:numzers+1,zeros(1,numzers+2),'-','Color',[.85 .85 .85],'LineWidth',0.5)
orders = params.aberrations(:,1:2);
allxticks = 1:numzers;
allxticklabels = cell(numzers,1);
for jzer = 1:numzers
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end
errorbar(1:numzers,mean(theta(6:end,:),2),std(theta(6:end,:),1,2));
xticks(allxticks)
xtickangle(25)
xticklabels(allxticklabels)
xlim([0 numzers+1])
ylim([-80 80])
yticks(-80:40:80)
yticklabels(-80:40:80)
xlabel('zernike mode (n,m)');
ylabel('rms value (m\lambda)');

% plot fitted zernike surfaces and data points (black dots)
figure
    set(gcf,'Position',[25 90 1255 880])
for nn = 1:size(zers,2)
    subplot(3,4,nn)
    surf(xim,yim,Zsurface(:,:,nn))
    hold on; axis square; shading interp;
    scatter3(xn,yn,zers(:,nn),2,'o','filled','MarkerFaceColor',[0 0 0])
    [r2,rmse] = rsquare(zers(:,nn),Zpredict(:,nn));
    hTitle=title(strcat(['Znm(' allxticklabels{nn} '), R2 = ' num2str(r2,2) ', RMSE = ' num2str(rmse,2)]));
    xlabel('x (\mu m)');
    ylabel('y (\mu m)');
    zlabel('rms value (m\lambda)');
    xlim([-60 60]);
    ylim([-60 60]);
    xticks(-60:30:60)
    yticks(-60:30:60)
    
end