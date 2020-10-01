function [thetastore,mustore,dmudthetastore,meritstore,numiters] = localization(allspots,theta0,params)
fprintf('\nStart fitting %i instances:\n',params.Ncfg); tic;
% This function finds the parameters for a single 2D-image assuming
% vectorial PSF models.
%

% parameter settings
Ncfg = params.Ncfg;
varfit = params.varfit;
tollim = params.tollim;
Nitermax = params.Nitermax;
numparams = params.numparams;

% pre-allocation
numiters = zeros(1,Ncfg);
mustore = zeros(params.Mx,params.My,params.K,Ncfg);
dmudthetastore = zeros(params.Mx,params.My,params.K,numparams,Ncfg);
thetastore = zeros(numparams,Ncfg,Nitermax+1);
meritstore = zeros(Ncfg,Nitermax+1);

flg_nat = params.flg_nat;
if flg_nat
    natPredictions = params.natPredictions;
else
    natPredictions = zeros(Ncfg,1);
end

% setup parallel loop
if params.flg_parallel
    p = gcp('nocreate');
    if isempty(p)
        parpool;
    end
    parforArg = Inf;
else
    parforArg = 0;
end

parfor (jcfg = 1:Ncfg, parforArg)
    % for jcfg = 1:Ncfg
    
    if flg_nat
        [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params,natPredictions(jcfg,:));
    else
        [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);
    end
    
    % pre-allocate
    thetatemp = zeros(numparams,Nitermax+1);
    merittemp = zeros(1,Nitermax+1);
    
    % initial values and max/min values
    theta = theta0(:,jcfg)';
    spots = allspots(:,:,:,jcfg);
    
    [thetamin,thetamax] = thetalimits(params,theta);
    
    thetaretry = (thetamax+thetamin)/2;
    
    [mu,dmudtheta] = poissonrate(params,theta,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
    [merit,grad,Hessian] = likelihood(params,spots,mu,dmudtheta,varfit);
    
    thetatemp(:,1) = theta;
    merittemp(1) = merit;
    meritprev = merit;
    
    % start iteration loop
    iiter = 1;
    monitor = 2*tollim;
    alamda = 1;
    alamdafac = 10;
    
    while ((iiter<=Nitermax) && (monitor>tollim))
        
        % check for det(H)=0 in order to avoid inversion of H
        matty = Hessian+alamda*diag(diag(Hessian));
        if (abs(det(matty))>2*eps)
            
            % update parameters
            thetatry = thetaupdate(theta,thetamax,thetamin,thetaretry,grad,Hessian,alamda,params);
            
            % calculate update merit function
            [mu,dmudtheta] = poissonrate(params,thetatry,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
            [merittry,gradtry,Hessiantry] = likelihood(params,spots,mu,dmudtheta,varfit);
            dmerit = merittry-merit;
            
            % modify Levenberg-Marquardt parameter
            if (dmerit<0)
                alamda = alamdafac*alamda;
            else
                alamda = alamda/alamdafac;
                theta = thetatry;
                merit = merittry;
                grad = gradtry;
                Hessian = Hessiantry;
                dmerit = merit-meritprev;
                monitor = abs(dmerit/merit);
                meritprev = merit;
                thetaretry = theta;
            end
            
        else
            alamda = alamdafac*alamda;
        end
        
        % store values and update counter
        thetatemp(:,iiter+1) = theta;
        merittemp(iiter+1) = merit;
        
        iiter = iiter+1; % update counter
        
    end
    
    % store values
    numiters(jcfg) = iiter-1;
    
    for jiter = iiter+1:Nitermax+1
        merittemp(jiter) = merit;
        thetatemp(:,jiter) = theta;
    end
    
    mustore(:,:,:,jcfg) = mu;
    dmudthetastore(:,:,:,:,jcfg) = dmudtheta;
    thetastore(:,jcfg,:) = thetatemp;
    meritstore(jcfg,:) = merittemp;
    
    % print update
    if rem(jcfg,round(Ncfg/10)) == 0
        fprintf('fitting instance # %i...\n',jcfg)
    end
    
end

% add offset to merit function
meritoffset = meritoffsetcalc(allspots,params.varfit);
for jiter = 1:Nitermax+1
    meritstore(:,jiter) = meritstore(:,jiter)+meritoffset;
end

% print run time
fprintf(['\nMLE fit routine (spot/second): ' num2str(toc,3) 's (' num2str(params.Ncfg/toc,5) ')\n'])

end

%%% tmp fun
% calculate theta estimation limits (max/min values)
function [thetamin,thetamax] = thetalimits(params,theta)

fitmodel = params.fitmodel;
zernikecoefsmax = 0.25*params.lambda*ones(1,size(params.aberrations,1));

roisizex = params.Mx*params.pixelsize;
roisizey = params.My*params.pixelsize;
xmin = -roisizex/2;
xmax = roisizex/2;
ymin = -roisizey/2;
ymax = roisizey/2;
zmin = params.zspread(1);
zmax = params.zspread(2);
azimmax = 2*pi;
polamax = pi;

if strcmp(fitmodel,'xy')
    thetamin = [xmin,ymin,theta(3)/10,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3)];
elseif strcmp(fitmodel,'xy-azim')
    thetamin = [xmin,ymin,theta(3)/10,0,0];
    thetamax = [xmax,ymax,roisizey/2,2*theta(3),theta(3),azimmax];
elseif strcmp(fitmodel,'xy-azim-pola')
    thetamin = [xmin,ymin,theta(3)/10,0,0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,polamax];
elseif strcmp(fitmodel,'xy-azim-pola-diffusion')
    thetamin = [xmin,ymin,theta(3)/10,0,0,-0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,polamax,1];
end

if strcmp(fitmodel,'xyz')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4)];
elseif strcmp(fitmodel,'xyz-azim-pola')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax];
elseif strcmp(fitmodel,'xy-azim-diffusion')
    thetamin = [xmin,ymin,theta(3)/10,0,0,0];
    thetamax = [xmax,ymax,2*theta(3),theta(3),azimmax,1];
elseif strcmp(fitmodel,'xyz-azim-pola-diffusion')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,0];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,1];
elseif strcmp(fitmodel,'xyz-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),zernikecoefsmax];
elseif strcmp(fitmodel,'xyz-azim-pola-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,zernikecoefsmax];
elseif strcmp(fitmodel,'xyz-azim-pola-diffusion-aberrations')
    thetamin = [xmin,ymin,zmin,theta(4)/10,0,0,0,0,-zernikecoefsmax];
    thetamax = [xmax,ymax,zmax,2*theta(4),theta(4),azimmax,polamax,1,zernikecoefsmax];
end

end