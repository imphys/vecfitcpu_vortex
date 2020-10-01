function [thetafinal,thetamean,thetastd,crlbmean] = get_statistics(params,object,theta,crlb,outliers)
% calculation of the mean and standard deviation of the estiamtes

% error
thetafinal = theta-object;

if strcmp(params.fitmodel,'xyz-azim-pola-diffusion')
    
    % convert orientations to half-sphere
    tmpcfg = theta(end-2,:)>pi;
    theta(end-2,tmpcfg) = theta(end-2,tmpcfg)-pi;
    theta(end-1,tmpcfg) = pi-theta(end-1,tmpcfg);
    
    % calculate orientation error using unit vectors
    Aest = Avec(theta(end-2,:),theta(end-1,:));
    Aobj = Avec(object(end-2,:),object(end-1,:));
    
    flip = sum(sign(Aobj)==sign(Aest))<2; % flipped match
    
    norm_azim = mod(theta(end-2,:)-object(end-2,:),pi);
    absdiff_azim = min(pi-norm_azim,norm_azim);
    absdiff_azim(norm_azim>pi/2) = -absdiff_azim(norm_azim>pi/2);
    
    norm_pola = mod(theta(end-1,:)-object(end-1,:),pi);
    norm_pola(flip) = mod(pi-theta(end-1,flip)-object(end-1,flip),pi);
    absdiff_pola = min(pi-norm_pola,norm_pola);
    absdiff_pola(norm_pola>pi/2) = -absdiff_pola(norm_pola>pi/2);
    
    thetafinal(end-2,:) = absdiff_azim;
    thetafinal(end-1,:) = absdiff_pola;
   
    outliers = outliers|(theta(7,:)>160*pi/180)|(theta(7,:)<20*pi/180);
    
elseif strcmp(params.fitmodel,'xy-azim-pola')
    
    % convert orientations to half-sphere
    tmpcfg = theta(end-1,:)>pi;
    theta(end-1,tmpcfg) = theta(end-1,tmpcfg)-pi;
    theta(end,tmpcfg) = pi-theta(end,tmpcfg);
    
    % calculate orientation error using unit vectors
    Aest = Avec(theta(end-1,:),theta(end,:));
    Aobj = Avec(object(end-1,:),object(end,:));
    
    flip = sum(sign(Aobj)==sign(Aest))<2; % flipped match
    
    norm_azim = mod(theta(5,:)-object(5,:),pi);
    absdiff_azim = min(pi-norm_azim,norm_azim);
    absdiff_azim(norm_azim>pi/2) = -absdiff_azim(norm_azim>pi/2);
    
    norm_pola = mod(theta(6,:)-object(6,:),pi);
    norm_pola(flip) = mod(pi-theta(6,flip)-object(6,flip),pi);
    absdiff_pola = min(pi-norm_pola,norm_pola);
    absdiff_pola(norm_pola>pi/2) = -absdiff_pola(norm_pola>pi/2);
    
    thetafinal(5,:) = absdiff_azim;
    thetafinal(6,:) = absdiff_pola;
    
    outliers = outliers|(theta(6,:)>160*pi/180)|(theta(6,:)<20*pi/180);
    
end

thetastd = std(thetafinal(:,~outliers),0,2);
thetamean = mean(thetafinal(:,~outliers),2);
crlbmean = mean(crlb(:,~outliers),2);

end

% calculate spherical vector (A)
function [A] = Avec(azim,pola)

A = [cos(azim).*sin(pola); sin(azim).*sin(pola); cos(pola)];

end


