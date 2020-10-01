function outliers = get_outliers(theta,meritstore,numiters,params)
% outlier removal based on statistics of maximum log-likelihood, on spot
% size, and on found position

Nitermax = params.Nitermax;

meritfinal = squeeze(meritstore(:,end))';
meanmeritfinal = mean(meritfinal);
stdmeritfinal = std(meritfinal);

outliers = (abs(theta(1,:))>3*params.pixelsize)|(abs(theta(2,:))>3*params.pixelsize);
outliers = outliers|(meritfinal<meanmeritfinal-3*stdmeritfinal);
outliers = outliers|(numiters==Nitermax);

fprintf(['\nNumber of outliers: ' num2str(sum(outliers))])
fprintf(['\nNumber of unconverged spots: ' num2str(sum(numiters==params.Nitermax)) '\n']);