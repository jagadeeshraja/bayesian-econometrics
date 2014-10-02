
%%% R Jagadeesh 10550
%load in the data set. Here use house price data from hprice.txt
load hprice.txt;
n=size(hprice,1);
y=hprice(:,1);
x=hprice(:,2:5);
x=[ones(n,1) x];
k=5;

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=5000;
b0(4,1)=10000;
b0(5,1)=10000;
s02=1/4.0e-8;
capv0=2.4*eye(k);
capv0(2,2)=6e-7;
capv0(3,3)=.15;
capv0(4,4)=.6;
capv0(5,5)=.6;
capv0inv=inv(capv0);

%Call script which carries actually does posterior analysis
post;

%save the log of marginal likelihood for later use
lmargun=lmarglik;

%Print out whatever you want
'Hyperparameters for informative natural conjugate prior'
b0
capv0
v0
s02

'Posterior results based on Informative Prior'
b1
bsd
probpos
bhpdi95
bhpdi99
hmean
hsd
lmarglik
ystarm
ystarsd
ystarcapv



%Hyperparameters for noninformative prior
v0=0;
capv0inv=0*eye(k);


%Call script which carries actually does posterior analysis
post;

%Print out whatever you want
'Posterior results based on Noninformative Prior'
b1
bsd
probpos
bhpdi95
bhpdi99
hmean
hsd
ystarm
ystarsd
ystarcapv

%posterior odds ratio
%evaluate log of marginal likelihood for restricted models with beta(j)=0
%this involves posterior inference for each of 5 restricted models
%there are better ways of programming this, but I simply do posterior
%analysis here for one model at a time
postodds=zeros(k,1);
x=hprice(:,2:5);
k=4;

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(1,1)=10;
b0(2,1)=5000;
b0(3,1)=10000;
b0(4,1)=10000;
s02=1/4.0e-8;
capv0=2.4*eye(k);
capv0(1,1)=6e-7;
capv0(2,2)=.15;
capv0(3,3)=.6;
capv0(4,4)=.6;
capv0inv=inv(capv0);

%Call script which carries actually does posterior analysis
post;

postodds(1,1)=exp(lmarglik-lmargun);

x=hprice(:,3:5);
x=[ones(n,1) x];
k=4;

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=5000;
b0(3,1)=10000;
b0(4,1)=10000;
s02=1/4.0e-8;
capv0=2.4*eye(k);
capv0(2,2)=.15;
capv0(3,3)=.6;
capv0(4,4)=.6;
capv0inv=inv(capv0);

%Call script which carries actually does posterior analysis
post;

postodds(2,1)=exp(lmarglik-lmargun);

x1=hprice(:,2);
x2=hprice(:,4:5);
x=[ones(n,1) x1 x2];
k=4;

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=10000;
b0(4,1)=10000;
s02=1/4.0e-8;
capv0=2.4*eye(k);
capv0(2,2)=6e-7;
capv0(3,3)=.6;
capv0(4,4)=.6;
capv0inv=inv(capv0);

%Call script which carries actually does posterior analysis
post;

postodds(3,1)=exp(lmarglik-lmargun);

x1=hprice(:,2:3);
x2=hprice(:,5);
x=[ones(n,1) x1 x2];
k=4;

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=5000;
b0(4,1)=10000;
s02=1/4.0e-8;
capv0=2.4*eye(k);
capv0(2,2)=6e-7;
capv0(3,3)=.15;
capv0(4,4)=.6;
capv0inv=inv(capv0);

%Call script which carries actually does posterior analysis
post;

postodds(4,1)=exp(lmarglik-lmargun);

x=hprice(:,2:4);
x=[ones(n,1) x];
k=4;

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=5000;
b0(4,1)=10000;
s02=1/4.0e-8;
capv0=2.4*eye(k);
capv0(2,2)=6e-7;
capv0(3,3)=.15;
capv0(4,4)=.6;
capv0inv=inv(capv0);

%Call script which carries actually does posterior analysis
post;

postodds(5,1)=exp(lmarglik-lmargun);

postodds

%posterior model probabilities for each restricted model
modelprob = postodds./(1+postodds)
