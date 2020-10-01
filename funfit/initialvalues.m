function [thetainit,errorfun] = initialvalues(allspots,params)
fprintf('\nInitial values (spots/second): '); tic;
% This function provides initial values for the fit parameters by means of
% a centroid fit.
%

K = params.K;
Mx = params.Mx;
My = params.My;
Ncfg = params.Ncfg;
zmin = params.zspread(1);
zmax = params.zspread(2);
deltaz = zmax-zmin;
fitmodel = params.fitmodel;
numparams = params.numparams;

if strcmp(params.doetype,'none')
    photonflux = 2.5; % rough photon flux correction
else
    photonflux = 2.5; % rough photon flux correction
end

% sampling coordinates in image plane
[~,~,XX,YY] = get_coords(params);

errorfun = zeros(Ncfg,1);
thetainit = zeros(numparams,Ncfg);
for jcfg = 1:Ncfg
    
    % rename config
    fxyk = allspots(:,:,:,jcfg);
    if strcmp(params.excitation,'zstack')
        fxy = fxyk(:,:,round(K/2));
    else
        fxy = sum(fxyk,3);
    end
    
    % background estimate from the median value of the rim pixels
    fxytmp = fxy';
    rimpx = [fxytmp(1:Mx) fxytmp(end-(Mx-1):end) fxy(1:My) fxy(end-(My-1):end)];
    bg = min(max(median(rimpx),1),50);
    % estimate of signal photon count
    
    if strcmp(params.excitation,'zstack')
        fxy = fxyk(:,:,round(K/2))-bg;
    else
        fxyk = max(fxyk-bg/K,1e-3);
        fxy = sum(fxyk,3);
    end
    
    % raw moments
    m00 = sum(fxy,1:2);
    m10 = sum(XX.*fxy,1:2);
    m01 = sum(YY.*fxy,1:2);
    m20 = sum(XX.^2.*fxy,1:2);
    m02 = sum(YY.^2.*fxy,1:2);
    m11 = sum(XX.*YY.*fxy,1:2);
    
    % centroids (lateral)
    xc = m10/m00;
    yc = m01/m00;
    
    % central moments
    mu20 = m20-m00*xc^2;
    mu11 = m11-m00*xc*yc;
    mu02 = m02-m00*yc^2;
    mu21 = m20-m20*yc-2*m11*yc+2*m00*xc^2*yc;
    mu12 = m02-m02*xc-2*m11*yc+2*m00*xc*yc^2;
    
    % centroid estimate of axial position
    if contains(fitmodel,'xyz')
        z0 = 1250*m11/(m20+m02);
        z0 = min(z0,zmax+0.5*deltaz);
        z0 = max(z0,zmin-0.5*deltaz);
        z0 = 0;
    end
    if contains(fitmodel,'pola')
        eig1 = ((mu20+mu02)+sqrt((mu20-mu02)^2+4*mu11^2))/2;
        eig2 = ((mu20+mu02)-sqrt((mu20-mu02)^2+4*mu11^2))/2;
        eccentricity = sqrt(1-eig2/eig1);
        pola0 = mod(pi*eccentricity,pi);
        % pola0 = pi/4;
    end
    
    if K == 1
        switch params.doetype
            case 'none'
                azim0 = pi/4;
            case 'vortex'
                azim0 = mod(-atan2(mu12,mu21),pi);
                % azim0 = pi/4;
        end
        
    elseif K>1
        azim0 = pi/4;
        pola0 = pi/4;
    end
    
    g20 = .75; % .499;
    
    switch fitmodel
        case 'xy'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = m00*photonflux;
            thetainit(4,jcfg) = bg;
        case 'xyz'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = m00*photonflux;
            thetainit(5,jcfg) = bg;
        case 'xyz-azim'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = m00*photonflux;
            thetainit(5,jcfg) = bg;
            thetainit(6,jcfg) = azim0;
        case 'xy-azim'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = m00*photonflux;
            thetainit(4,jcfg) = bg;
            thetainit(5,jcfg) = azim0;
        case 'xy-azim-diffusion'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = m00*photonflux;
            thetainit(4,jcfg) = bg;
            thetainit(5,jcfg) = azim0;
            thetainit(6,jcfg) = g20;
        case 'xyz-azim-pola'
            thetainit(1,jcfg) = yc;
            thetainit(2,jcfg) = xc;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = m00*photonflux;
            thetainit(5,jcfg) = bg;
            thetainit(6,jcfg) = azim0;
            thetainit(7,jcfg) = pola0;
        case 'xy-azim-pola'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = m00*photonflux;
            thetainit(4,jcfg) = bg;
            thetainit(5,jcfg) = azim0;
            thetainit(6,jcfg) = pola0;
        case 'xy-azim-pola-diffusion'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = m00*photonflux;
            thetainit(4,jcfg) = bg;
            thetainit(5,jcfg) = azim0;
            thetainit(6,jcfg) = pola0;
            thetainit(7,jcfg) = g20;
        case 'xyz-azim-pola-diffusion'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = m00*photonflux;
            thetainit(5,jcfg) = bg;
            thetainit(6,jcfg) = azim0;
            thetainit(7,jcfg) = pola0;
            thetainit(8,jcfg) = g20;
            
        case 'xyz-aberrations'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = m00*photonflux;
            thetainit(5,jcfg) = bg;
            thetainit(6:params.numparams,jcfg) = 0;
            
        case 'xyz-azim-pola-aberrations'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = m00*photonflux;
            thetainit(5,jcfg) = bg;
            thetainit(6,jcfg) = azim0;
            thetainit(7,jcfg) = pola0;
            thetainit(8:params.numparams,jcfg) = 0;
            
        case 'xyz-azim-pola-diffusion-aberrations'
            thetainit(1,jcfg) = xc;
            thetainit(2,jcfg) = yc;
            thetainit(3,jcfg) = z0;
            thetainit(4,jcfg) = m00*photonflux;
            thetainit(5,jcfg) = bg;
            thetainit(6,jcfg) = azim0;
            thetainit(7,jcfg) = pola0;
            thetainit(8,jcfg) = g20;
            thetainit(9:params.numparams,jcfg) = 0;
            
    end
end

fprintf([num2str(toc,3) 's (' num2str(params.Ncfg/toc,5) ')\n'])

