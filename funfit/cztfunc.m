function dataout = cztfunc(datain,Amt,Bmt,Dmt,params)
% This function evaluates the FT via the czt-algorithm
% arguments: datain = input data, dimensions K x N
%            A,B,D = auxiliary vectors computed in prechirpz, must have
%            lengths N, M, and L=N+M-1
% function value: dataout = output data, dimensions K x M

N = params.cztN; % size(Amt,2);
M = params.cztM; % size(Bmt,2);
L = params.cztL; % N+M-1;

cztin = zeros(size(datain,1),L);
% cztin(:,1:N)= Amt.*datain;
cztin(:,1:N)= datain;
temp = Dmt.*fft(cztin,[],2);
cztout = ifft(temp,[],2);
dataout = Bmt.*cztout(:,1:M);
  


