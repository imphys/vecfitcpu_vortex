function [logL,gradlogL,HessianlogL] = likelihood(params,image,mu,dmudtheta,varfit)
% returns the log-likelihood, as well as first and second order
% derivatives w.r.t. the parameters for a noisy image M, measured in
% number of photons per pixel, and Poisson-rate mu with first order
% derivatives dmudtheta. The log-likelihood is Poisson+readout-noise based.
%

numparams = params.numparams;

% calculation of weight factors
keps = 1e3*eps;
mupos = double(mu>0).*mu + double(mu<0)*keps;
weight = (image-mupos)./(mupos+varfit);
dweight = (image+varfit)./(mupos+varfit).^2;

% log-likelihood, gradient vector and Hessian matrix
logL = sum((image+varfit).*log(mupos+varfit)-(mupos+varfit),1:3);
gradlogL = permute(sum(weight.*dmudtheta,1:3),[5 4 3 2 1]);
HessianlogL = zeros(numparams,numparams);
for ii = 1:numparams
    for jj = 1:numparams
        HessianlogL(ii,jj) = sum(-dweight.*dmudtheta(:,:,:,ii).*dmudtheta(:,:,:,jj),1:3);
    end
end





