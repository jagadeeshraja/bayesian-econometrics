
function [z] = rnddtruncnorm(mu,sigma,a,b)         % rnd_dtrunc_norm

% Generates random variables from a Doubly Truncated Normal Distribution
% mu    : mean
% sigma : standard deviation
% a     : left truncated point
% b     : right truncated point

u = rand;
w = normcdf((a - mu)/sigma) + u*(normcdf((b - mu)/sigma) - normcdf((a - mu)/sigma));
z = mu + sigma*norminv(w);