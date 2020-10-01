function meritoffset = meritoffsetcalc(allspots,varfit)
% This function calculates the merit function offset for an accurate
% determination of the log-likelihood.
%

[Mx,My,K,Ncfg] = size(allspots);
meritoffset = zeros(Ncfg,1,'Like',allspots);
for jcfg = 1:Ncfg
    dummat = allspots(:,:,:,jcfg);
    % set negative and zero pixels values to one/10 to avoid log-singularity
    dummat = max(dummat,ones(size(dummat))/10);
    meritoffset(jcfg) = 0;
    for ii=1:Mx
        for jj = 1:My
            for kk = 1:K
                meritoffset(jcfg) = meritoffset(jcfg)-gammln(dummat(ii,jj,kk)+1+varfit);
            end
        end
    end
end

