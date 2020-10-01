function thetanew = thetaupdate(thetaold,thetamax,thetamin,thetaretry,grad,Hessian,alamda,params)
% This function calculates the new parameters in each iteration step.
%

fitmodel = params.fitmodel;

% update of fit parameters via Levenberg-Marquardt
Bmat = Hessian+alamda*diag(diag(Hessian));
dtheta = -Bmat\grad';
thetanew = thetaold+dtheta';

% enforce physical boundaries in angular space
if contains(fitmodel,'azim-pola')
    
    if contains(fitmodel,'xyz')
        if thetanew(6)>thetamax(6)
            thetanew(6) = thetanew(6)-pi;
            thetanew(7) = pi-thetanew(7);
        elseif thetanew(6)<thetamin(6)
            thetanew(6) = thetanew(6)+pi;
            thetanew(7) = pi-thetanew(7);
        end
        
        thetanew(6) = mod(thetanew(6),2*pi);
        thetanew(7) = mod(thetanew(7),pi);
        
        if contains(fitmodel,'diffusion')
            if thetanew(8)>thetamax(8)
                thetanew(8) = thetamax(8)-thetanew(8);
            elseif thetanew(8)<thetamin(8)
                thetanew(8) = thetamin(8)-thetanew(8);
            end
        end
        
    elseif contains(fitmodel,'xy')
        
        if thetanew(5)>thetamax(5)
            thetanew(5) = thetanew(5)-pi;
            thetanew(6) = pi-thetanew(6);
        elseif thetanew(5)<thetamin(5)
            thetanew(5) = thetanew(5)+pi;
            thetanew(6) = pi-thetanew(6);
        end
        
        thetanew(5) = mod(thetanew(5),2*pi);
        thetanew(6) = mod(thetanew(6),pi);
        
    end
    
end

% enforce physical boundaries in parameter space.
for jj=1:length(thetaold)
    if ((thetanew(jj)>thetamax(jj))||(thetanew(jj)<thetamin(jj)))
        thetanew(jj) = thetaretry(jj);
    end
end
