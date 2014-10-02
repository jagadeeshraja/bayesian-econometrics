
function [z] = drawlatent(y,x,beta)

% The function draws latent variables in Bayesian estimation of a
% bivariate probit model

n = size(y,1);
z = zeros(n,1);

for i=1:n
    if y(i)==1
       z(i,1)=rnddtruncnorm(x(i,:)*beta,1,0,Inf);
    else
       z(i,1)=rnddtruncnorm(x(i,:)*beta,1,-Inf,0);
    end
end