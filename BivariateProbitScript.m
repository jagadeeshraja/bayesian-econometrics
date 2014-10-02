%%% R Jagadeesh 10550
% Bayesian Estimation of a Bivariate Probit Model

data    = xlsread('binary.csv');

% naming the variables

admit   = data(:,1);
gre     = data(:,2);
gpa     = data(:,3);
rank    = data(:,4);

rank2   = logical(rank==2);
rank3   = logical(rank==3);
rank4   = logical(rank==4);

n       = size(admit,1);                            % no. of observations
y       = admit;
x       = [ones(n,1) gre gpa rank2 rank3 rank4];
k       = size(x,2);                                % no. of covariates

%% Estimation

% Prior Specification
b0      = zeros(k,1);
B0      = eye(k);
invB0   = eye(k);                % inv(B0) = B0

% Starting values
beta    = zeros(k,1);

% Posterior Variance
Btilde  = inv(B0 + x'*x);        % inv(B0) = B0, posterior variance should 
                                 % be calculated outside the loop to speed the algorithm 
nsim    = 10000;
burn    = 0.25*nsim;

storebetas = zeros(k,nsim);
% Gibbs sampling
tic
h = waitbar(0,'Simulation in Progress');
for i = 1:nsim
    % Draw z from the full conditional posterior
    z = drawlatent(y,x,beta);
    % Draw beta from the full conditional posterior
    beta = drawbeta(z,x,b0,invB0,Btilde);
    
    storebetas(:,i) = beta;
    waitbar(i/nsim);
end
close(h)

% Solution
postmeanbeta = mean(storebetas(:, burn+1:nsim),2)
%    -1.6754  0.0012  0.2928    -0.4163   -0.7837   -0.9231

poststdbeta  = std(storebetas(:, burn+1:nsim)')
%    0.5572    0.0006    0.1728    0.1857    0.1993    0.2346







