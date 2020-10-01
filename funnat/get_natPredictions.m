function [prediction] = get_natPredictions(xn,yn,RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6)

% Z(5,6)
prediction(:,1) = 2*RAstig3(1)*(xn.^3.*yn+xn.*yn.^3) + RAstig3(2)*(xn.^2.*yn+yn.^3) +...
    RAstig3(3)*(xn.^3+xn.*yn.^2) + RAstig3(4)*(xn.^2+yn.^2) +...
    RAstig3(6)*(xn.*yn.^2-xn.^3) + RAstig3(7)*(2*xn.*yn.^2) +...
    RAstig3(8)*yn - RAstig3(9).*xn + RAstig3(10)*(2*xn.*yn) +...
    RAstig3(11)*xn+RAstig3(12)*yn+RAstig3(13);

prediction(:,2) = RAstig3(1)*(yn.^4-xn.^4)-RAstig3(2)*(xn.^3+xn.*yn.^2)+ ...
    RAstig3(3)*(xn.^2.*yn+yn.^3)+RAstig3(5)*(xn.^2+yn.^2)-...
    RAstig3(6)*(2*xn.^2.*yn)+RAstig3(7)*(yn.^3-xn.^2.*yn)+RAstig3(8).*xn+...
    RAstig3(9)*yn+RAstig3(10)*(yn.^2-xn.^2)+RAstig3(11)*yn-...
    RAstig3(12)*xn+RAstig3(14);

% Z(7,8)
prediction(:,3) = RComa3(1)*(xn.^3+xn.*yn.^2)+RComa3(2)*xn.^2+RComa3(3)*xn.*yn+...
    RComa3(4)*xn+RComa3(5)*(xn.^2+yn.^2)+RComa3(7)*yn+RComa3(8)*xn+RComa3(9);
prediction(:,4) = RComa3(1)*(yn.^3+xn.^2.*yn)+RComa3(2)*xn.*yn+RComa3(3)*yn.^2+...
    RComa3(4)*yn+RComa3(6)*(xn.^2+yn.^2)+RComa3(7)*xn-RComa3(8)*yn+RComa3(10);

% Z(9,10)
prediction(:,5) = RTrefoil(1)*(3*yn.^2.*xn-xn.^3)+RTrefoil(2)*(yn.^2-xn.^2)+...
    RTrefoil(3)*(2*xn.*yn)+RTrefoil(4)*xn+RTrefoil(5)*yn+RTrefoil(6);
prediction(:,6) = RTrefoil(1)*(yn.^3-3*xn.^2.*yn)-RTrefoil(2)*(2*xn.*yn)+...
    RTrefoil(3)*(yn.^2-xn.^2)+RTrefoil(4)*yn-RTrefoil(5)*xn+RTrefoil(7);

% Z11surf
prediction(:,7) = RCurv5(1)*(xn.^2+yn.^2)+RCurv5(2)*xn+RCurv5(3)*yn+RCurv5(4);

% Z(12,13)
prediction(:,8) = RAstig5(1)*(2*xn.*yn)+RAstig5(2)*yn+RAstig5(3)*xn+RAstig5(4);
prediction(:,9) = RAstig5(1)*(yn.^2-xn.^2)-RAstig5(2)*xn+RAstig5(3)*yn+RAstig5(5);

% Z(16,17)
prediction(:,10) = RComa5(1)*xn+RComa5(2);
prediction(:,11) = RComa5(1)*yn+RComa5(3);

% Z22surf
prediction(:,12) = RCurv6(1)*(xn.^2+yn.^2)+RCurv6(2)*xn+RCurv6(3)*yn+RCurv6(4);

% Z22surf
% prediction(:,13) = RCurv8(1)*(xn.^2+yn.^2)+RCurv8(2)*xn+RCurv8(3)*yn+RCurv8(4);

% Z22surf
% prediction(:,14) = RCurv10(1)*(xn.^2+yn.^2)+RCurv10(2)*xn+RCurv10(3)*yn+RCurv10(4);