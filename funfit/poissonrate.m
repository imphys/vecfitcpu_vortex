function [mu,dmudtheta] = poissonrate(params,theta,PupilMatrix,allzernikes,wavevector,wavevectorzimm)
% Returns the Poisson-rates for all pixels and all first order derivatives
% w.r.t. the parameters theta.

K = params.K;
% m = params.m;
Mx = params.Mx;
My = params.My;

fitmodel = params.fitmodel;
numparams = params.numparams;
% alpha = params.alpha;
% beta = params.beta;

switch fitmodel
    case 'xy'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        
    case 'xyz'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        
    case 'xy-azim'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        azim = theta(5);
        params.azim = azim;
        
    case 'xy-azim-pola'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        azim = theta(5);
        pola = theta(6);
        params.azim = azim;
        params.pola = pola;
        
    case 'xyz-azim-pola'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        
    case 'xy-azim-diffusion'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        azim = theta(5);
        params.azim = azim;
        params.g2 = theta(6);
        
    case 'xy-azim-pola-diffusion'
        params.xemit = theta(1);
        params.yemit = theta(2);
        Nph = theta(3);
        Nbg = theta(4);
        azim = theta(5);
        pola = theta(6);
        params.azim = azim;
        params.pola = pola;
        params.g2 = theta(7);
        
    case 'xyz-azim-pola-diffusion'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        params.g2 = theta(8);
        
    case 'xyz-aberrations'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        params.aberrations(:,3) = theta(6:end);
        
    case 'xyz-azim-pola-aberrations'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        params.aberrations(:,3) = theta(8:end);
        
    case 'xyz-azim-pola-diffusion-aberrations'
        params.xemit = theta(1);
        params.yemit = theta(2);
        params.zemit = theta(3);
        Nph = theta(4);
        Nbg = theta(5);
        azim = theta(6);
        pola = theta(7);
        params.azim = azim;
        params.pola = pola;
        params.g2 = theta(8);
        params.aberrations(:,3) = theta(9:end);
        
end

switch params.excitation
    case 'constant'
        P = 1;
        dPdazim = 0;
        dPdpola = 0;
    case 'zstack'
        P = ones(1,K);
        dPdazim = zeros(1,K);
        dPdpola = zeros(1,K);
end

% update pupil function
if contains(fitmodel,'aberrations')
    [wavevector,wavevectorzimm,~,allzernikes,PupilMatrix] = get_pupil_matrix(params);
end

[FieldMatrix,FieldMatrixDerivatives] = get_field_matrix_derivatives(params,PupilMatrix,allzernikes,wavevector,wavevectorzimm);
[PSF,PSFder] = get_psfs_derivatives(params,PupilMatrix,FieldMatrix,FieldMatrixDerivatives);

% get Poisson rate and derivatives
mu = zeros(Mx,My,K);
dmudtheta = zeros(Mx,My,K,numparams);
if strcmp(params.excitation,'zstack')
    
    if strcmp(fitmodel,'xyz-aberrations')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1:3) = Nph*PSFder(:,:,:,1:3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1;
        dmudtheta(:,:,:,6:end) = Nph*PSFder(:,:,:,4:end);
    end
    
    if strcmp(fitmodel,'xyz-azim-pola-aberrations')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1) = Nph*PSFder(:,:,:,1);
        dmudtheta(:,:,:,2) = Nph*PSFder(:,:,:,2);
        dmudtheta(:,:,:,3) = Nph*PSFder(:,:,:,3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1/K;
        dmudtheta(:,:,:,6) = Nph*PSFder(:,:,:,4);
        dmudtheta(:,:,:,7) = Nph*PSFder(:,:,:,5);
        dmudtheta(:,:,:,8:end) = Nph*PSFder(:,:,:,6:end);
    end
    
    if strcmp(fitmodel,'xyz-azim-pola-diffusion-aberrations')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1) = Nph*PSFder(:,:,:,1);
        dmudtheta(:,:,:,2) = Nph*PSFder(:,:,:,2);
        dmudtheta(:,:,:,3) = Nph*PSFder(:,:,:,3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1/K;
        dmudtheta(:,:,:,6) = Nph*PSFder(:,:,:,4);
        dmudtheta(:,:,:,7) = Nph*PSFder(:,:,:,5);
        dmudtheta(:,:,:,8) = Nph*PSFder(:,:,:,6);
        dmudtheta(:,:,:,9:end) = Nph*PSFder(:,:,:,7:end);
    end
    
    if strcmp(fitmodel,'xyz-azim-pola-diffusion')
        mu(:,:,:) = Nph*PSF(:,:,:)+Nbg;
        dmudtheta(:,:,:,1) = Nph*PSFder(:,:,:,1);
        dmudtheta(:,:,:,2) = Nph*PSFder(:,:,:,2);
        dmudtheta(:,:,:,3) = Nph*PSFder(:,:,:,3);
        dmudtheta(:,:,:,4) = PSF;
        dmudtheta(:,:,:,5) = 1/K;
        dmudtheta(:,:,:,6) = Nph*PSFder(:,:,:,4);
        dmudtheta(:,:,:,7) = Nph*PSFder(:,:,:,5);
        dmudtheta(:,:,:,8) = Nph*PSFder(:,:,:,6);
    end
    
else
    for k = 1:K
        % PSF model
        mu(:,:,k) = Nph*P(k)*PSF+Nbg/K;
        
        % get derivatives of Poisson rate w.r.t. fit parameters
        if strcmp(fitmodel,'xy-azim-pola') && K>1
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
            dmudtheta(:,:,k,5) = Nph*P(k)*PSFder(:,:,3)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdpola(k)*PSF;
        end
        
        if strcmp(fitmodel,'xy')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
        end
        
        if strcmp(fitmodel,'xyz')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = Nph*P(k)*PSFder(:,:,3);
            dmudtheta(:,:,k,4) = P(k)*PSF;
            dmudtheta(:,:,k,5) = 1/K;
        end
        
        if strcmp(fitmodel,'xy-azim')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
            dmudtheta(:,:,k,5) = Nph*dPdazim(k)*PSF+Nph*P(k)*PSFder(:,:,3);
        end
        
        if strcmp(fitmodel,'xyz-azim-pola')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = Nph*P(k)*PSFder(:,:,3);
            dmudtheta(:,:,k,4) = P(k)*PSF;
            dmudtheta(:,:,k,5) = 1/K;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,7) = Nph*P(k)*PSFder(:,:,5)+Nph*dPdpola(k)*PSF;
        end
        
        if strcmp(fitmodel,'xy-azim-diffusion')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
            dmudtheta(:,:,k,5) = Nph*P(k)*PSFder(:,:,3)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4);
        end
        
        if strcmp(fitmodel,'xy-azim-pola-diffusion')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = P(k)*PSF;
            dmudtheta(:,:,k,4) = 1/K;
            dmudtheta(:,:,k,5) = Nph*P(k)*PSFder(:,:,3)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdpola(k)*PSF;
            dmudtheta(:,:,k,7) = Nph*P(k)*PSFder(:,:,5);
        end
        
        if strcmp(fitmodel,'xyz-azim-pola-diffusion')
            dmudtheta(:,:,k,1) = Nph*P(k)*PSFder(:,:,1);
            dmudtheta(:,:,k,2) = Nph*P(k)*PSFder(:,:,2);
            dmudtheta(:,:,k,3) = Nph*P(k)*PSFder(:,:,3);
            dmudtheta(:,:,k,4) = P(k)*PSF;
            dmudtheta(:,:,k,5) = 1/K;
            dmudtheta(:,:,k,6) = Nph*P(k)*PSFder(:,:,4)+Nph*dPdazim(k)*PSF;
            dmudtheta(:,:,k,7) = Nph*P(k)*PSFder(:,:,5)+Nph*dPdpola(k)*PSF;
            dmudtheta(:,:,k,8) = Nph*P(k)*PSFder(:,:,6);
        end
        
    end
end