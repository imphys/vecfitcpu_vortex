function [merit,grad,Hessian] = meritfun(theta,spotscfg,parameters)
% This function computes the merit function for the fmincon implementation
% of the vectorial MLE routine. It calculates the Poisson-rates and
% derivatives, and subsequently the log-likelihood and derivatives.
%

[mu,dmudtheta] = poissonrate(theta,parameters);
[merit,grad,Hessian] = likelihood(spotscfg,mu,dmudtheta,parameters.varfit);
merit = -merit;
grad = -grad;
Hessian = -Hessian;

