function [FieldMatrix,FieldMatrixDerivatives] = get_field_matrix_derivatives(params,PupilMatrix,allzernikes,wavevector,wavevectorzimm)
% This function calculates the field matrix A_{jk}, which gives the j-th
% electric field component proportional to the k-th dipole vector
% component, as well as the derivatives of A_{jk} w.r.t. the xyz coordinates
% of the emitter and w.r.t. the emission wavelength lambda.
%
% parameters: NA, refractive indices of medium, wavelength (in nm),
% nominal emitter position (in nm) with z-position from
% cover slip-medium interface, spot footprint (in nm), axial range (in nm),
% sampling in pupil with (even), sampling in image plane (odd), sampling in
% axial direction

K = params.K;
Mx = params.Mx;
My = params.My;
% Mz = params.Mz;
xemit = params.xemit;
yemit = params.yemit;
zemit = params.zemit;
Npupil = params.Npupil;
fitmodel = params.fitmodel;

% wavevector = params.wavevector;
% PupilMatrix = params.PupilMatrix;
% wavevectorzimm = params.wavevectorzimm;

% calculate auxiliary vectors for chirpz
Ax = params.Axmt;
Bx = params.Bxmt;
Dx = params.Dxmt;
Ay = params.Aymt;
By = params.Bymt;
Dy = params.Dymt;

% calculation Zernike mode normalization
if contains(fitmodel,'aberrations')
    %         allzernikes = params.allzernikes;
    orders = params.aberrations(:,1:2);
    normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
end

% preallocate
if contains(fitmodel,'xyz')
    numders = 3;
elseif contains(fitmodel,'xy')
    numders = 2;
end

if strcmp(params.excitation,'zstack')
    
    % define z-stack
    zmin = params.zrange(1);
    zmax = params.zrange(2);
    ImageSizez = (zmax-zmin)/2;
    DzImage = 2*ImageSizez/K;
    ZImage = zmin+DzImage/2:DzImage:zmax;
    %     DzImage = 2*ImageSizez/(K-1);
    %     ZImage = (zmin:DzImage:zmax);
    
    % preallocate
    FieldMatrix = zeros(Mx,My,K,2,3);
    FieldMatrixDerivatives = zeros(Mx,My,K,numders,2,3);
    PupilFunctionDerivatives = zeros(Npupil,Npupil,numders,2,3);
    
    % run over z-stack
    for jz = 1:numel(ZImage)
        zemitrun = ZImage(jz);
        
        % phase contribution due to position of the emitter
        Wlateral = xemit*wavevector(:,:,1)+yemit*wavevector(:,:,2);
        if strcmp(params.ztype,'stage')
            Wpos = Wlateral+(zemit-zemitrun)*wavevectorzimm;
        elseif strcmp(params.ztype,'medium')
            Wpos = Wlateral+(zemit-zemitrun)*wavevector(:,:,3);
        end
        PositionPhaseMask = exp(-1i*Wpos);
        
        % Pupil function
        PupilFunction = PositionPhaseMask.*PupilMatrix;
        
        % pupil functions for xy-derivatives
        PupilFunctionDerivatives(:,:,1,:,:) = -1i*wavevector(:,:,1).*PupilFunction;
        PupilFunctionDerivatives(:,:,2,:,:) = -1i*wavevector(:,:,2).*PupilFunction;
        
        % pupil functions for z-derivatives
        if contains(fitmodel,'xyz')
            if strcmp(params.ztype,'stage')
                PupilFunctionDerivatives(:,:,3,:,:) = -1i*wavevectorzimm.*PupilFunction;
            elseif strcmp(params.ztype,'medium')
                PupilFunctionDerivatives(:,:,3,:,:) = -1i*wavevector(:,:,3).*PupilFunction;
            end
        end
        
        % Field matrix and derivatives
        IntermediateImage = cztfunc2D(PupilFunction,Ay,By,Dy,params);
        FieldMatrix(:,:,jz,:,:) = cztfunc2D(IntermediateImage,Ax,Bx,Dx,params);
        IntermediateImage = cztfunc3D(PupilFunctionDerivatives,Ay,By,Dy,params);
        FieldMatrixDerivatives(:,:,jz,1:numders,:,:) = cztfunc3D(IntermediateImage,Ax,Bx,Dx,params);
        
        % pupil functions for Zernike mode-derivative and FT to matrix elements
        if contains(fitmodel,'aberrations')
            for jzer = 1:size(params.aberrations,1)
                jder = numders+jzer;
                PupilFunction = (2*pi*1i*normfac(jzer)*allzernikes(:,:,jzer)/params.lambda).*PositionPhaseMask.*PupilMatrix;
                IntermediateImage = cztfunc2D(PupilFunction,Ay,By,Dy,params);
                FieldMatrixDerivatives(:,:,jz,jder,:,:) = cztfunc2D(IntermediateImage,Ax,Bx,Dx,params);
            end
        end
        
    end
    
else
    
    % phase contribution due to position of the emitter
    Wlateral = xemit*wavevector(:,:,1)+yemit*wavevector(:,:,2);
    if strcmp(params.ztype,'stage')
        Wpos = Wlateral+zemit*wavevectorzimm;
    elseif strcmp(params.ztype,'medium')
        Wpos = Wlateral+zemit*wavevector(:,:,3);
    end
    PositionPhaseMask = exp(-1i*Wpos);
    
    % Pupil function and derivatives
    PupilFunction = PositionPhaseMask.*PupilMatrix;
    
    % preallocate
    PupilFunctionDerivatives = zeros(Npupil,Npupil,numders,2,3);
    % pupil functions for xy-derivatives
    PupilFunctionDerivatives(:,:,1,:,:) = -1i*wavevector(:,:,1).*PupilFunction;
    PupilFunctionDerivatives(:,:,2,:,:) = -1i*wavevector(:,:,2).*PupilFunction;
    % pupil functions for z-derivatives
    if contains(fitmodel,'xyz')
        if strcmp(params.ztype,'stage')
            PupilFunctionDerivatives(:,:,3,:,:) = -1i*wavevectorzimm.*PupilFunction;
        else
            PupilFunctionDerivatives(:,:,3,:,:) = -1i*wavevector(:,:,3).*PupilFunction;
        end
    end
    
    % Field matrix and derivatives
    IntermediateImage = cztfunc2D(PupilFunction,Ay,By,Dy,params);
    FieldMatrix = cztfunc2D(IntermediateImage,Ax,Bx,Dx,params);
    IntermediateImage = cztfunc3D(PupilFunctionDerivatives,Ay,By,Dy,params);
    FieldMatrixDerivatives = cztfunc3D(IntermediateImage,Ax,Bx,Dx,params);
    
end

%% plotting intermediate results
if params.debugmode
    jz = ceil(Mz/2);
    figure
    for itel = 1:2
        for jtel = 1:3
            tempim = FieldMatrix(:,:,itel,jtel,jz);
            subplot(2,3,3*(itel-1)+jtel)
            imagesc(abs(tempim))
            title(strcat('amplitude i=',num2str(itel),', j=',num2str(jtel)))
            axis square
            axis off
        end
    end
    figure
    for itel = 1:2
        for jtel = 1:3
            tempim = FieldMatrix(:,:,itel,jtel,jz);
            subplot(2,3,3*(itel-1)+jtel)
            imagesc(angle(tempim)*180/pi)
            title(strcat('phase i=',num2str(itel),', j=',num2str(jtel)))
            axis square
            axis off
        end
    end
end
