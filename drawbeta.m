
function [beta] = drawbeta(z,x,b0,invB0,Btilde)

% The function draws beta in Bayesian estimation of a
% bivariate probit model

betatilde   = Btilde*(invB0*b0 + x'*z);

beta        = mvnrnd(betatilde,Btilde)';