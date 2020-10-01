function theta = roi2fov(theta,roixy,params)

% roi to fov coordinates
xtmp = theta(1,:)/params.pixelsize;
ytmp = theta(2,:)/params.pixelsize;
theta(1,:) = roixy(:,1)'+ytmp;
theta(2,:) = roixy(:,2)'+xtmp;

% sphere2hemisphere  (0<azim<pi, 0<pola<pi)
if contains(params.fitmodel,'xyz')
    if contains(params.fitmodel,'azim-pola')
        tmpcfg = theta(6,:)>pi;
        theta(6,tmpcfg) = theta(6,tmpcfg)-pi;
        theta(7,tmpcfg) = pi-theta(7,tmpcfg);
    end
else
    if contains(params.fitmodel,'azim-pola')
        tmpcfg = theta(5,:)>pi;
        theta(5,tmpcfg) = theta(5,tmpcfg)-pi;
        theta(6,tmpcfg) = pi-theta(6,tmpcfg);
    end
end