function [RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6] = get_natCoefficients(xn,yn,zers,params)

Ncfg = params.Ncfg;
tol = 1e-6;

xn2 = xn.^2;
xn3 = xn.^3;
xn4 = xn.^4;
yn2 = yn.^2;
yn3 = yn.^3;
yn4 = yn.^4;

% B = ones(Ncfg,1);
Z5 = zers(:,1);
Z6 = zers(:,2);
Z7 = zers(:,3);
Z8 = zers(:,4);
Z9 = zers(:,5);
Z10 = zers(:,6);
Z11 = zers(:,7);
Z12 = zers(:,8);
Z13 = zers(:,9);
Z16 = zers(:,10);
Z17 = zers(:,11);
Z22 = zers(:,12);

MAstig3 = zeros(2*Ncfg,14);
MComa3 = zeros(2*Ncfg,10);
MTrefoil = zeros(2*Ncfg,7);
MCurv5 = zeros(Ncfg,4);
MCurv6 = zeros(Ncfg,4);
MAstig5 = zeros(2*Ncfg,5);
MComa5 = zeros(2*Ncfg,3);
for nn = 1:Ncfg
    mm = 2*nn;
    
    x = [xn(nn) xn2(nn) xn3(nn) xn4(nn)];
    y = [yn(nn) yn2(nn) yn3(nn) yn4(nn)];
    
    Z56(mm-1:mm,1) = [Z5(nn); Z6(nn)];
    Z78(mm-1:mm,1) = [Z7(nn); Z8(nn)];
    Z910(mm-1:mm,1) = [Z9(nn); Z10(nn)];
    Z1213(mm-1:mm,1) = [Z12(nn);Z13(nn)];
    Z1617(mm-1:mm,1) = [Z16(nn);Z17(nn)];
    
    %LS Matrix for Astig3 (third order astigmatism)
    M1 = [2*(x(3)*y(1)+x(1)*y(3)) x(2)*y(1)+y(3) x(3)+x(1)*y(2) x(2)+y(2) 0 ...
        x(1)*y(2)-x(3) 2*x(1)*y(2) y(1) -x(1) 2*x(1)*y(1) ...
        x(1) y(1) 1 0];
    M2 = [y(4)-x(4) -x(3)-x(1)*y(2) x(2)*y(1)+y(3) 0 x(2)+y(2) ...
        -2*x(2)*y(1) y(3)-x(2)*y(1) x(1) y(1) y(2)-x(2) ...
        y(1) -x(1) 0 1];
    MAstig3(mm-1:mm,:) = [M1;M2];
    
    %LS Matrix for Coma3 (field cubed coma: W331)
    M1 = [x(3)+x(1)*y(2) x(2) x(1)*y(1) x(1) x(2)+y(2) 0 y(1) x(1) 1 0];
    M2 = [y(3)+x(2)*y(1) x(1)*y(1) y(2) y(1) 0 x(2)+y(2) x(1) -y(1) 0 1];
    MComa3(mm-1:mm,:) = [M1;M2];
    
    %LS Matrix for Trefoil5
    M1 = [3*y(2)*x(1)-x(3) y(2)-x(2) 2*x(1)*y(1) x(1) y(1) 1 0];
    M2 = [y(3)-3*x(2)*y(1) -2*x(1)*y(1) y(2)-x(2) y(1) -x(1) 0 1];
    MTrefoil(mm-1:mm,:) = [M1;M2];
    
    % (Z40) LS Matrix for Curv5 (fifth order field curvature: W240)
    MCurv5(nn,:) = [x(2)+y(2) x(1) y(1) 1];
    
    % (Z60) LS Matrix for 5th order spherical
    MCurv6(nn,:) = [x(2)+y(2) x(1) y(1) 1];
    
    %LS Matrix for Astig3 (third order astigmatism)
    M1 = [2*x(1)*y(1) y(1) x(1) 1 0];
    M2 = [y(2)-x(2) -x(1) y(1) 0 1];
    MAstig5(mm-1:mm,:) = [M1;M2];
    
    %LS Matrix for Coma5 (fifth order coma: W151)
    M1 = [x(1) 1 0];
    M2 = [y(1) 0 1];
    MComa5(mm-1:mm,:) = [M1;M2];
    
end

RAstig3 = pinv(MAstig3,tol)*Z56;
fAstig3 = MAstig3*RAstig3;

RComa3 = pinv(MComa3,tol)*Z78;
fComa3 = MComa3*RComa3;

RTrefoil = pinv(MTrefoil,tol)*Z910;
fTrefoil = MTrefoil*RTrefoil;

RCurv5 = pinv(MCurv5,tol)*Z11;
fCurv5 = MCurv5*RCurv5;

RAstig5 = pinv(MAstig5,tol)*Z1213;
fAstig5 = MAstig5*RAstig5;

RComa5 = pinv(MComa5,tol)*Z1617;
fComa5 = MComa5*RComa5;

RCurv6 = pinv(MCurv6,tol)*Z22;
fCurv6 = MCurv6*RCurv6;

% R^2 values
[r25,rmse5] = rsquare(Z5,fAstig3(1:2:end));
[r26,rmse6] = rsquare(Z6,fAstig3(2:2:end));
[r27,rmse7] = rsquare(Z7,fComa3(1:2:end));
[r28,rmse8] = rsquare(Z8,fComa3(2:2:end));
[r29,rmse9] = rsquare(Z9,fTrefoil(1:2:end));
[r210,rmse10] = rsquare(Z10,fTrefoil(2:2:end));
[r211,rmse11] = rsquare(Z11,fCurv5);
[r212,rmse12] = rsquare(Z12,fAstig5(1:2:end));
[r213,rmse13] = rsquare(Z13,fAstig5(2:2:end));
[r216,rmse16] = rsquare(Z16,fComa5(1:2:end));
[r217,rmse17] = rsquare(Z17,fComa5(2:2:end));
[r222,rmse22] = rsquare(Z22,fCurv6);

% display R^2
disp(' ')
disp(['Astig3:  R2/RMSE = ' num2str(r25,2) '/' num2str(rmse5,2) ', ' num2str(r26,2) '/' num2str(rmse6,2)])
disp(['Coma3:   R2/RMSE = ' num2str(r27,2) '/' num2str(rmse7,2) ', ' num2str(r28,2) '/' num2str(rmse8,2)])
disp(['Trefoil: R2/RMSE = ' num2str(r29,2) '/' num2str(rmse9,2) ', ' num2str(r210,2) '/' num2str(rmse10,2)])
disp(['Curv5:   R2/RMSE = ' num2str(r211,2) '/' num2str(rmse11,2)])
disp(['Astig5:  R2/RMSE = ' num2str(r212,2) '/' num2str(rmse12,2) ', ' num2str(r213,2) '/' num2str(rmse13,2)])
disp(['Coma5:   R2/RMSE = ' num2str(r216,2) '/' num2str(rmse16,2) ', ' num2str(r217,3) '/' num2str(rmse17,2)])
disp(['Curv6:   R2/RMSE = ' num2str(r222,2) '/' num2str(rmse22,2)])
