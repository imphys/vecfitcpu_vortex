function [normint_free,normint_fixed] = get_normalization(params,PupilMatrix)
% This function computes the PSF normalization by evaluating the energy
% flow through the lens aperture
%

% parameters
NA = params.NA;
pola = params.pola;
azim = params.azim;
lambda = params.lambda;
Npupil = params.Npupil;
pixelsize = params.pixelsize;

% dipole unit vector
dipor = [sin(pola)*cos(azim) sin(pola)*sin(azim) cos(pola)];

% intensity matrix
IntensityMatrix = zeros(3,3);
for itel = 1:3
    for jtel = 1:3
        %                 for ztel = 1:2
        %                     pupmat1 = PupilMatrix(:,:,ztel,itel);
        %                     pupmat2 = PupilMatrix(:,:,ztel,jtel);
        %                     IntensityMatrix(itel,jtel) = IntensityMatrix(itel,jtel)+sum(sum(real(pupmat1.*conj(pupmat2))));
        %                 end
        
        pupmat1 = PupilMatrix(:,:,:,itel);
        pupmat2 = PupilMatrix(:,:,:,jtel);
        IntensityMatrix(itel,jtel) = sum(real(pupmat1.*conj(pupmat2)),'all');
    end
end

% normalization to take into account discretization correctly
DxyPupil = 2*NA/lambda/Npupil;
normfac = DxyPupil^2/(pixelsize)^2;
IntensityMatrix = normfac*IntensityMatrix;

% evaluation normalization factors
normint_free = sum(diag(IntensityMatrix))/3;
normint_fixed = dipor*(IntensityMatrix*transpose(dipor));

% ring average for emission profile as a function of radial pupil coordinate
if params.debugmode
    offset = [0,0];
    pixelsizes = [1 1]*2/Npupil;
    numbins = Npupil;
    radial_intensity_matrix = zeros(3,3,numbins);
    for itel = 1:3
        for jtel = 1:3
            for ztel = 1:2
                pupmat1 = PupilMatrix(:,:,ztel,itel);
                pupmat2 = PupilMatrix(:,:,ztel,jtel);
                Amatin = real(pupmat1.*conj(pupmat2));
                [~,Avecout,~,radiusvec] = radialavgmat(Amatin,numbins,offset,pixelsizes);
                radial_intensity_matrix(itel,jtel,:) = squeeze(radial_intensity_matrix(itel,jtel,:))+Avecout;
            end
        end
    end
    radial_intensity_free = zeros(numbins,1);
    radial_intensity_fixed = zeros(numbins,1);
    for jbin = 1:numbins
        tempmat = squeeze(radial_intensity_matrix(:,:,jbin));
        radial_intensity_free(jbin) = sum(diag(tempmat))/3;
        radial_intensity_fixed(jbin) = dipor*(tempmat*transpose(dipor));
    end
    
    figure
    hold on
    box on
    plot(radiusvec,radial_intensity_free,'r')
    plot(radiusvec,radial_intensity_fixed,'b')
    legend('free','fixed','Location','NorthWest')
    xlim([0 1])
    xlabel('normalized pupil coordinate')
    ylabel('emitted intensity (a.u.)')
    
    figure
    hold on
    box on
    angle_imm = asin(radiusvec*params.NA/params.refimm)';
    plot(angle_imm*180/pi,radial_intensity_free.*cos(angle_imm),'r')
    plot(angle_imm*180/pi,radial_intensity_fixed.*cos(angle_imm),'b')
    legend('free','fixed','Location','NorthWest')
    xlim([0 90])
    xlabel('Angle in immersion fluid (deg)')
    ylabel('emitted intensity (a.u.)')
    
end