function errM = get_fiterror(mu,spots,params)
% fprintf('\nchi-square: '); tic;
% This function calculates the normalized chisquare value
% fprintf([num2str(toc,3) 's\n']);

normalize = params.K*params.Mx*params.My;
mu = double(mu>0).*mu + double(mu<0)*1e3*eps;

if strcmp(params.excitation,'zstack') || params.K>1
    
    mu = permute(squeeze(mu),[4 1 2 3]);
    spots = permute(squeeze(spots),[4 1 2 3]);
    
    chi2 = sum((mu-spots).^2./mu ,2:4);
    LLR = 2*sum(mu - spots + spots.*log(spots) - spots.*log(mu) ,2:4);
    SSE = sum(mu-spots,2:4).^2';
    
else

    mu = permute(squeeze(mu),[3 1 2]);
    spots = permute(squeeze(spots),[3 1 2]);
    
    chi2 = sum((mu-spots).^2./mu ,2:3);
    LLR = 2*sum(mu - spots + spots.*log(spots) - spots.*log(mu) ,2:3);
    SSE = sum(mu-spots,2:3).^2;
    
end

errM(1,:) = chi2;
errM(2,:) = LLR;
errM(3,:) = SSE;
errM = errM/normalize;