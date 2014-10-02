%%% R Jagadeesh
%For a given data set and values for prior hyperparameters
%Calculate OLS and posterior quantities, marginal likelihood
%and predictive (Predictive sets x-star=.5

%Ordinary least squares quantities
bols = inv(x'*x)*x'*y;
s2 = (y-x*bols)'*(y-x*bols)/(n-k);
bolscov = s2*inv(x'*x);
bolssd=zeros(k,1);
for i = 1:k
bolssd(i,1)=sqrt(bolscov(i,i));
end
v=n-k;

%Posterior hyperparameters for Normal-Gamma
xsquare=x'*x;
v1=v0+n;
capv1inv = capv0inv+ xsquare;
capv1=inv(capv1inv);
b1 = capv1*(capv0inv*b0 + xsquare*bols);
if det(capv0inv)>0
    v1s12 = v0*s02 + v*s2 + (bols-b0)'*inv(capv0 + inv(xsquare))*(bols-b0);
else
    v1s12 = v0*s02 + v*s2;
end
s12 = v1s12/v1;

bcov = capv1*v1s12/(v1-2);
bsd=zeros(k,1);
for i = 1:k
bsd(i,1)=sqrt(bcov(i,i));
end

%posterior probability that each element of beta is positive
%as well as 95 and 99 HPDIs for each element of beta

probpos=zeros(k,1);
bhpdi95=zeros(k,2);
bhpdi99=zeros(k,2);

%get quantiles of t for calculating HPDIs
invcdf95=tdis_inv(.975,v1);
invcdf99=tdis_inv(.995,v1);

for i = 1:k
    tnorm = -b1(i,1)/sqrt(s12*capv1(i,i));
    probpos(i,1) = 1 - tdis_cdf(tnorm,v1);
    bhpdi95(i,1) = b1(i,1)-invcdf95*sqrt(s12*capv1(i,i));
    bhpdi95(i,2) = b1(i,1)+invcdf95*sqrt(s12*capv1(i,i));
    bhpdi99(i,1) = b1(i,1)-invcdf99*sqrt(s12*capv1(i,i));
    bhpdi99(i,2) = b1(i,1)+invcdf99*sqrt(s12*capv1(i,i));
end



%posterior mean and variance of error precision
hmean = 1/s12;
hvar=2/(v1s12);
hsd=sqrt(hvar);

%predictive inference
if k == 5
xstar = [1 5000 2 2 1];
ystarm = xstar*b1;
ystarcapv = (1+ xstar*capv1*xstar')*s12;
ystarv = ystarcapv*v1/(v1-2);
ystarsd=sqrt(ystarv);
end

%log of marginal likelihood for the model if prior is informative
if det(capv0inv)>0;
    intcon=gammaln(.5*v1) + .5*v0*log(v0*s02)- gammaln(.5*v0) -.5*n*log(pi);
    lmarglik=intcon + .5*log(det(capv1)/det(capv0)) - .5*v1*log(v1s12);
end