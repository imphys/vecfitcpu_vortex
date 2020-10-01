function [wavevector,wavevectorzimm,Waberration,allzernikes,PupilMatrix] = get_pupil_matrix(params,natPredictions)
% This function calculates the pupil matrix Q_{jk}, which gives the j-th
% electric field component proportional to the k-th dipole vector
% component.
%

NA = params.NA;
refmed = params.refmed;
refcov = params.refcov;
refimm = params.refimm;
refimmnom = params.refimmnom;
lambda = params.lambda;
Npupil = params.Npupil;
aberrations = params.aberrations;
zvals = params.zvals;

if nargin>1
    aberrations(:,3) = natPredictions';
end

% pupil radius (in diffraction units) and pupil coordinate sampling
PupilSize = 1.0;
DxyPupil = 2*PupilSize/Npupil;
XYPupil = -PupilSize+DxyPupil/2:DxyPupil:PupilSize;
[YPupil,XPupil] = meshgrid(XYPupil,XYPupil);

% calculation of relevant Fresnel-coefficients for the interfaces
% between the medium and the cover slip and between the cover slip
% and the immersion fluid
% The Fresnel-coefficients should be divided by the wavevector z-component
% of the incident medium, this factor originates from the
% Weyl-representation of the emitted vector spherical wave of the dipole.
% CosThetaMed = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);

% SAF correction 26-07-2020
% Clockwise rotation (sqrt assumes counterclockwise; Eulers formula)
argMed = 1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2;
phiMed = atan2(0,argMed);
CosThetaMed = sqrt(abs(argMed)).*(cos(phiMed/2)-1j*sin(phiMed/2));

CosThetaCov = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refcov^2);
CosThetaImm = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refimm^2);
CosThetaImmnom = sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refimmnom^2);

FresnelPmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaCov+refcov*CosThetaMed);
FresnelSmedcov = 2*refmed*CosThetaMed./(refmed*CosThetaMed+refcov*CosThetaCov);
FresnelPcovimm = 2*refcov*CosThetaCov./(refcov*CosThetaImm+refimm*CosThetaCov);
FresnelScovimm = 2*refcov*CosThetaCov./(refcov*CosThetaCov+refimm*CosThetaImm);
FresnelP = FresnelPmedcov.*FresnelPcovimm;
FresnelS = FresnelSmedcov.*FresnelScovimm;

% setting of vectorial functions
Phi = atan2(YPupil,XPupil);
CosPhi = cos(Phi);
SinPhi = sin(Phi);
CosTheta = CosThetaMed; % sqrt(1-(XPupil.^2+YPupil.^2)*NA^2/refmed^2);
SinTheta = sqrt(1-CosTheta.^2);

pvec(:,:,1) = FresnelP.*CosTheta.*CosPhi;
pvec(:,:,2) = FresnelP.*CosTheta.*SinPhi;
pvec(:,:,3) = -FresnelP.*SinTheta;
svec(:,:,1) = -FresnelS.*SinPhi;
svec(:,:,2) = FresnelS.*CosPhi;
svec(:,:,3) = 0;

% PolarizationVector = zeros(Npupil,Npupil,2,3);
PolarizationVector(:,:,1,:) = CosPhi.*pvec-SinPhi.*svec;
PolarizationVector(:,:,2,:) = SinPhi.*pvec+CosPhi.*svec;

% definition aperture
ApertureMask = double((XPupil.^2+YPupil.^2)<1.0);

% aplanatic amplitude factor
% combining this factor with the Fresnel-coefficient factors T_{p} and T_{s}
% the amplitude scales as [sqrt(cos(theta_imm))/cos(theta_med)] x T_{p,s}
% where T_{p} = T_{p,med->cov} x T_{p,cov->imm}
% and T_{s} = T_{s,med->cov} x T_{s,cov->imm}
% with T_{p,med->cov} = 2*ref_med*cos(theta_med)/[ref_med*cos(theta_cov)+ref_cov*cos(theta_med)]
% and T_{s,med->cov} = 2*ref_med*cos(theta_med)/[ref_med*cos(theta_med)+ref_cov*cos(theta_cov)]
% in case of index matching the overall amplitude scaling is with
% [1/sqrt(cos(theta_med))] x T_{p,s}
Amplitude = ApertureMask.*sqrt(CosThetaImm)./(refmed*CosThetaMed);

% calculation aberration function
Waberration = zeros(size(XPupil));
orders = aberrations(:,1:2);
zernikecoefs = squeeze(aberrations(:,3));
normfac = sqrt(2*(orders(:,1)+1)./(1+double(orders(:,2)==0)));
zernikecoefs = normfac.*zernikecoefs;
allzernikes = get_zernikefunctions(orders,XPupil,YPupil);
for j = 1:numel(zernikecoefs)
    Waberration = Waberration+zernikecoefs(j)*allzernikes(:,:,j);
end

% calculation DOE/SLM phase function
if ~strcmp(params.doetype,'none')
    switch params.doetype
        case 'vortex'
            numzones = length(params.ringradius)-1;
            
            rho = sqrt(XPupil.^2+YPupil.^2);
            phi = atan2(YPupil,XPupil);
            
            zoneindex = zeros(size(rho));
            for jj = 1:numzones
                zoneindex = zoneindex+jj*double(rho>params.ringradius(jj)).*double(rho<=params.ringradius(jj+1));
            end
            
            if params.doelevels==0
                ZoneFunction = (2*zoneindex-1).*phi/(2*pi)+0.5;
            else
                ZoneFunction = round((params.doelevels)*((2*zoneindex-1).*phi/(2*pi)+0.5))/(params.doelevels);
            end
            
            if isfield(params,'doevortexflip')
                if params.doevortexflip==1
                    ZoneFunction=max(max(ZoneFunction))-ZoneFunction;
                end
            end
            
    end
    
    switch params.doetype
        case 'vortex'
            DOEaberration = params.doephasedepth*(ZoneFunction-floor(ZoneFunction));
    end
    
    Waberration = Waberration+DOEaberration;
    
end

Waberration = Waberration+zvals(1)*refimm*CosThetaImm-zvals(2)*refimmnom*CosThetaImmnom-zvals(3)*refmed*CosThetaMed;
PhaseFactor = exp(2*pi*1i*Waberration/lambda);
Waberration = Waberration.*ApertureMask;

% compute pupil matrix
PupilMatrix = Amplitude.*PhaseFactor.*PolarizationVector;

% calculate wavevector inside immersion fluid and z-component inside medium
wavevector(:,:,1) = (2*pi*NA/lambda)*XPupil;
wavevector(:,:,2) = (2*pi*NA/lambda)*YPupil;
wavevector(:,:,3) = (2*pi*refmed/lambda)*CosThetaMed;
wavevectorzimm = (2*pi*refimm/lambda)*CosThetaImm;

% plotting intermediate results
if params.debugmode
    figure
    for itel = 1:2
        for jtel = 1:3
            tempim = PupilMatrix(:,:,itel,jtel);
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
            tempim = PupilMatrix(:,:,itel,jtel);
            subplot(2,3,3*(itel-1)+jtel)
            imagesc(angle(tempim)*180/pi)
            title(strcat('phase i=',num2str(itel),', j=',num2str(jtel)))
            axis square
            axis off
        end
    end
end

end