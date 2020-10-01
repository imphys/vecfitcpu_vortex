function dataout = cztfunc3D(datain,Amt,Bmt,Dmt,params)
% This function evaluates the FT via the czt-algorithm
% arguments: datain = input data, dimensions K x N
%            A,B,D = auxiliary vectors computed in prechirpz, must have
%            lengths N, M, and L=N+M-1
% function value: dataout = output data, dimensions K x M


N = params.cztN;
M = params.cztM;
L = params.cztL;

if contains(params.fitmodel,'xyz')
    dim = 3;
elseif contains(params.fitmodel,'xy')
    dim = 2;
end
    
cztin = zeros(size(datain,1),L,dim,2,3);

cztin(:,1:N,:,:,:) = Amt.*datain;
temp = Dmt.*fft(cztin,[],2);
cztout = ifft(temp,[],2);
dataout = Bmt.*cztout(:,1:M,:,:,:);
dataout = permute(dataout,[2 1 3 4 5]);