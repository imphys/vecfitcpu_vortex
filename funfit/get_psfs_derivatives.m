function [PSF,PSFder] = get_psfs_derivatives(params,PupilMatrix,FieldMatrix,FieldMatrixDerivatives)
% This function calculates the free or fixed dipole PSFs given the field
% matrix, the dipole orientation, and the pupil polarization, as well as
% the derivatives w.r.t. the xyz coordinates of the emitter and w.r.t. the
% emission wavelength lambda.
%
% parameters: emitter/absorber dipole orientation (characterized by angles
% pola and azim).
%

K = params.K;
Mx = params.Mx;
My = params.My;
g2 = params.g2;
azim = params.azim;
pola = params.pola;
fitmodel = params.fitmodel;

% PSF normalization
[normint_free,normint_fixed] = get_normalization(params,PupilMatrix);

% dipole vector and derivatives
dipor = [sin(pola)*cos(azim) sin(pola)*sin(azim) cos(pola)];
if contains(fitmodel,'azim')
    diporDerivatives(1,:) = [-sin(pola)*sin(azim) sin(pola)*cos(azim) 0];
end
if contains(fitmodel,'pola')
    diporDerivatives(2,:) = [cos(pola)*cos(azim) cos(pola)*sin(azim) -sin(pola)];
end

% parameter indices
if contains(fitmodel,'azim-pola')
    if contains(fitmodel,'xyz')
        nazim = 4;
        npola = 5;
    else
        nazim = 3;
        npola = 4;
    end
end

% calculation of free PSF and derivatives
if strcmp(params.dipoletype,'free') || strcmp(params.dipoletype,'diffusion')
    if strcmp(params.excitation,'zstack')
        
        FreePSF = 1/3*sum(abs(FieldMatrix).^2,4:5);
        tmpFieldMatrixDerivatives = permute(FieldMatrixDerivatives,[1:3 5:6 4]);
        tmpPSFderivatives = (2/3)*sum(real(conj(FieldMatrix).*tmpFieldMatrixDerivatives),4:5);
        FreePSFder = permute(tmpPSFderivatives,[1:3 6 4:5]);
        
        if strcmp(fitmodel,'xyz-azim-pola-aberrations')
            tmpFreePSFder = zeros(params.Mx,params.My,K,params.numparams-2);
            tmpFreePSFder(:,:,:,1:3) = FreePSFder(:,:,:,1:3);
            tmpFreePSFder(:,:,:,6:end) = FreePSFder(:,:,:,4:end);
            FreePSFder = tmpFreePSFder;
        end
        
        if strcmp(fitmodel,'xyz-azim-pola-diffusion-aberrations')
            tmpFreePSFder = zeros(params.Mx,params.My,K,params.numparams-2);
            tmpFreePSFder(:,:,:,1:3) = FreePSFder(:,:,:,1:3);
            tmpFreePSFder(:,:,:,7:end) = FreePSFder(:,:,:,4:end);
            FreePSFder = tmpFreePSFder;
        end
        
    else % ~zstack
        
        FreePSF = 1/3*sum(abs(FieldMatrix).^2,3:4);
        tmpFieldMatrixDerivatives = permute(FieldMatrixDerivatives,[1:2 4:5 3]);
        tmpPSFderivatives = (2/3)*sum(real(conj(FieldMatrix).*tmpFieldMatrixDerivatives),3:4);
        FreePSFder = permute(tmpPSFderivatives,[1:2 5 3:4]);
        
    end
    
    FreePSF = FreePSF/normint_free;
    FreePSFder = FreePSFder/normint_free;
    
    if strcmp(params.dipoletype,'free')
        PSF = FreePSF;
        PSFder = FreePSFder;
    end
    
end

% calculated of fixed PSF and derivatives
if strcmp(params.dipoletype,'fixed') || strcmp(params.dipoletype,'diffusion')
    if strcmp(params.excitation,'zstack')
        
        % electric field components and derivatives
        Ex = dipor(1)*FieldMatrix(:,:,:,1,1)+dipor(2)*FieldMatrix(:,:,:,1,2)+dipor(3)*FieldMatrix(:,:,:,1,3);
        Ey = dipor(1)*FieldMatrix(:,:,:,2,1)+dipor(2)*FieldMatrix(:,:,:,2,2)+dipor(3)*FieldMatrix(:,:,:,2,3);
        
        Exder = dipor(1)*FieldMatrixDerivatives(:,:,:,:,1,1)+dipor(2)*FieldMatrixDerivatives(:,:,:,:,1,2)+dipor(3)*FieldMatrixDerivatives(:,:,:,:,1,3);
        Eyder = dipor(1)*FieldMatrixDerivatives(:,:,:,:,2,1)+dipor(2)*FieldMatrixDerivatives(:,:,:,:,2,2)+dipor(3)*FieldMatrixDerivatives(:,:,:,:,2,3);
        
        if strcmp(fitmodel,'xyz-azim-pola-aberrations')
            tmpExder = zeros(params.Mx,params.My,K,params.numparams-2);
            tmpEyder = zeros(params.Mx,params.My,K,params.numparams-2);
            
            tmpExder(:,:,:,1:3) = Exder(:,:,:,1:3);
            tmpEyder(:,:,:,1:3) = Eyder(:,:,:,1:3);
            
            tmpExder(:,:,:,6:end) = Exder(:,:,:,4:end);
            tmpEyder(:,:,:,6:end) = Eyder(:,:,:,4:end);
            
            Exder = tmpExder;
            Eyder = tmpEyder;
        end
        
        if strcmp(fitmodel,'xyz-azim-pola-diffusion-aberrations')
            tmpExder = zeros(params.Mx,params.My,K,params.numparams-2);
            tmpEyder = zeros(params.Mx,params.My,K,params.numparams-2);
            
            tmpExder(:,:,:,1:3) = Exder(:,:,:,1:3);
            tmpEyder(:,:,:,1:3) = Eyder(:,:,:,1:3);
            
            tmpExder(:,:,:,7:end) = Exder(:,:,:,4:end);
            tmpEyder(:,:,:,7:end) = Eyder(:,:,:,4:end);
            
            Exder = tmpExder;
            Eyder = tmpEyder;
        end
        
        
        % azim
        if contains(fitmodel,'azim')
            Exder(:,:,:,nazim) = diporDerivatives(1,1)*FieldMatrix(:,:,:,1,1)+diporDerivatives(1,2)*FieldMatrix(:,:,:,1,2)+diporDerivatives(1,3)*FieldMatrix(:,:,:,1,3);
            Eyder(:,:,:,nazim) = diporDerivatives(1,1)*FieldMatrix(:,:,:,2,1)+diporDerivatives(1,2)*FieldMatrix(:,:,:,2,2)+diporDerivatives(1,3)*FieldMatrix(:,:,:,2,3);
        end
        
        % pola
        if contains(fitmodel,'pola')
            Exder(:,:,:,npola) = diporDerivatives(2,1)*FieldMatrix(:,:,:,1,1)+diporDerivatives(2,2)*FieldMatrix(:,:,:,1,2)+diporDerivatives(2,3)*FieldMatrix(:,:,:,1,3);
            Eyder(:,:,:,npola) = diporDerivatives(2,1)*FieldMatrix(:,:,:,2,1)+diporDerivatives(2,2)*FieldMatrix(:,:,:,2,2)+diporDerivatives(2,3)*FieldMatrix(:,:,:,2,3);
        end
        
        FixedPSF = abs(Ex).^2+abs(Ey).^2;
        FixedPSFder = 2*real(conj(Ex).*Exder)+2*real(conj(Ey).*Eyder);
        
    else
        
        Ex = dipor(1)*FieldMatrix(:,:,1,1)+dipor(2)*FieldMatrix(:,:,1,2)+dipor(3)*FieldMatrix(:,:,1,3);
        Ey = dipor(1)*FieldMatrix(:,:,2,1)+dipor(2)*FieldMatrix(:,:,2,2)+dipor(3)*FieldMatrix(:,:,2,3);
        Exder = dipor(1)*FieldMatrixDerivatives(:,:,:,1,1)+dipor(2)*FieldMatrixDerivatives(:,:,:,1,2)+dipor(3)*FieldMatrixDerivatives(:,:,:,1,3);
        Eyder = dipor(1)*FieldMatrixDerivatives(:,:,:,2,1)+dipor(2)*FieldMatrixDerivatives(:,:,:,2,2)+dipor(3)*FieldMatrixDerivatives(:,:,:,2,3);
        
        % azim
        if contains(fitmodel,'azim')
            Exder(:,:,nazim) = diporDerivatives(1,1)*FieldMatrix(:,:,1,1)+diporDerivatives(1,2)*FieldMatrix(:,:,1,2)+diporDerivatives(1,3)*FieldMatrix(:,:,1,3);
            Eyder(:,:,nazim) = diporDerivatives(1,1)*FieldMatrix(:,:,2,1)+diporDerivatives(1,2)*FieldMatrix(:,:,2,2)+diporDerivatives(1,3)*FieldMatrix(:,:,2,3);
            
        end
        % pola
        if contains(fitmodel,'pola')
            Exder(:,:,npola) = diporDerivatives(2,1)*FieldMatrix(:,:,1,1)+diporDerivatives(2,2)*FieldMatrix(:,:,1,2)+diporDerivatives(2,3)*FieldMatrix(:,:,1,3);
            Eyder(:,:,npola) = diporDerivatives(2,1)*FieldMatrix(:,:,2,1)+diporDerivatives(2,2)*FieldMatrix(:,:,2,2)+diporDerivatives(2,3)*FieldMatrix(:,:,2,3);
        end
        
        FixedPSF = abs(Ex).^2+abs(Ey).^2;
        FixedPSFder = 2*real(conj(Ex).*Exder)+2*real(conj(Ey).*Eyder);
        
    end
    
    FixedPSF = FixedPSF/normint_fixed;
    FixedPSFder = FixedPSFder/normint_fixed;
    
    if strcmp(params.dipoletype,'fixed')
        PSF = FixedPSF;
        PSFder = FixedPSFder;
    end
    
end

% calculated of diffusion PSF and derivatives
if strcmp(params.dipoletype,'diffusion')
    
    if strcmp(params.excitation,'zstack')
        if contains(fitmodel,'azim')
            FreePSFder(:,:,:,nazim) = zeros(Mx,My,K);
        end
        if contains(fitmodel,'pola')
            FreePSFder(:,:,:,npola) = zeros(Mx,My,K);
        end
    else
        if contains(fitmodel,'azim')
            FreePSFder(:,:,nazim) = zeros(Mx,My);
        end
        if contains(fitmodel,'pola')
            FreePSFder(:,:,npola) = zeros(Mx,My);
        end
    end
    
    PSF = (1-g2)*FreePSF+g2*FixedPSF;
    PSFder = (1-g2)*FreePSFder+g2*FixedPSFder;
    
    if strcmp(params.excitation,'zstack')
        
        if strcmp(fitmodel,'xyz-azim-pola-diffusion')
            PSFder(:,:,:,6) = -FreePSF+FixedPSF;
        elseif strcmp(fitmodel,'xy-azim-pola-diffusion')
            PSFder(:,:,:,5) = -FreePSF+FixedPSF;
        elseif strcmp(fitmodel,'xy-azim-diffusion')
            PSFder(:,:,:,4) = -FreePSF+FixedPSF;
        elseif strcmp(fitmodel,'xyz-azim-pola-diffusion-aberrations')
            PSFder(:,:,:,6) = -FreePSF+FixedPSF;
        end
        
    else
        
        if strcmp(fitmodel,'xyz-azim-pola-diffusion')
            PSFder(:,:,6) = -FreePSF+FixedPSF;
        elseif strcmp(fitmodel,'xy-azim-pola-diffusion')
            PSFder(:,:,5) = -FreePSF+FixedPSF;
        elseif strcmp(fitmodel,'xy-azim-diffusion')
            PSFder(:,:,4) = -FreePSF+FixedPSF;
        end
        
    end
    
end

% Blurs PSF due to the effect of a non-zero pixel size
% PSF = do_pixel_blurring(PSF,params);

% 3D convolution of the PSFs and derivatives with a bead
if isfield(params,'bead')
    if params.bead == true
        bead = create3DBead(params);
        PSF = convn(bead,PSF,'same');
        tempderivs = zeros(size(PSFder));
        if strcmp(params.excitation,'zstack')
            for jder = 1:size(PSFder,4)
                % tempderivs(:,:,:,jder) = convn(bead,squeeze(PSFder(:,:,:,jder)),'same');
                tempderivs(:,:,:,jder) = convn(bead,PSFder(:,:,:,jder),'same');
            end
        else
            for jder = 1:size(PSFder,3)
                % tempderivs(:,:,jder) = convn(bead,squeeze(PSFder(:,:,jder)),'same');
                tempderivs(:,:,jder) = convn(bead,PSFder(:,:,jder),'same');
            end
        end
        PSFder = tempderivs;
    end
end

end