function [A,B,D] = prechirpz(xsize,qsize,N,M)
% This function evaluates the auxiliary vectors for the evaluation of the
% FT via the czt-algorithm
% arguments: xsize = window in real space abs(x)<xsize
%            qsize = window in Fourier space abs(q)<qsize
%            N = # sample points in real space (even)
%            M = # sample points in Fourier space (odd)
% function value: A,B,D = auxiliary vectors of lengths N, M, and L=N+M-1

L = N+M-1;
sigma = -2*pi*xsize*qsize/N/M;
Gfac = (2*xsize/N)*exp(1i*sigma*(1-N)*(1-M));

% % computation of vectors via recursion
% Afac = exp(2*1i*sigma*(1-M));
% Bfac = exp(2*1i*sigma*(1-N));
% sqW = exp(2*1i*sigma);
% W = sqW^2;
% 
% Utmp = zeros(1,N);
% A = zeros(1,N);
% Utmp(1) = sqW*Afac;
% A(1) = 1.0;
% for i=2:N
%   A(i) = Utmp(i-1)*A(i-1);
%   Utmp(i) = Utmp(i-1)*W;
% end
% 
% Vtmp = zeros(1,M);
% B = ones(1,M);
% Vtmp(1) = sqW*Bfac;
% B(1) = Gfac;
% for i=2:M
%   B(i) = Vtmp(i-1)*B(i-1);
%   Vtmp(i) = Vtmp(i-1)*W;
% end
% 
% Utmp = zeros(1,max(N,M)+1);
% Vtmp = zeros(1,max(N,M)+1);
% Utmp(1) = sqW;
% Vtmp(1) = 1.0;
% for i=2:max(N,M)+1
%   Vtmp(i) = Utmp(i-1)*Vtmp(i-1);
%   Utmp(i) = Utmp(i-1)*W;
% end

% direct computation of vectors
allns = 0:(N-1);
A = exp(2*1i*sigma*allns.*(allns+1-M)); 
allms = 0:(M-1);
B = Gfac*exp(2*1i*sigma*allms.*(allms+1-N));
allnms = 0:max(N,M)+1;
Vtmp = exp(2*1i*sigma*allnms.^2);

D = ones(1,L);
for i=1:M
  D(i) = conj(Vtmp(i));
end
for i=1:N-1
  D(L+1-i) = conj(Vtmp(i+1));
end

D = fft(D);

end

