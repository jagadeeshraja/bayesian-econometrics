(Importance sampling)
%Initial Sample size
S=100;
%Posterior Distribution(Mean=0 And Variance=1)
f=@(x)(1/sqrt(2*pi)*exp(-x.^2/2));
%Proposal Distribution t3 Distribution
g=@(x)(6*sqrt(3))./(pi*(3+x.^2).^2);
% weight Fucntion 
w=@(x)f(x)./g(x);
%sampling value from proposal distribtuion.
theta_array = trnd(3,[S,1]);

%Weight Vector
weight=w(theta_array);
%Log of Wieght Vector;

log_w=log(weight);
% histogram plot
nbins = 25;
subplot(2,2,1);
hist(log_w,nbins)
%Plot for threshold values.
log_w(log_w<-0.19)=[];
subplot(2,2,2);
hist(log_w,nbins)
% Expectation
Expect_1= sum(theta_array.*weight)/sum(weight)

% Varaince
Vari_1= sum((theta_array.^2).*weight)/sum(weight) - (Expect_1^2)

%Initial Sample size
S=1000;
%Posterior Distribution(Mean=0 And Variance=1)
f=@(x)(1/sqrt(2*pi)*exp(-x.^2/2));
%Proposal Distribution t3 Distribution
g=@(x)(6*sqrt(3))./(pi*(3+x.^2).^2);
% weight Fucntion 
w=@(x)f(x)./g(x);
%sampling value from proposal distribtuion.
theta_array = trnd(3,[S,1]);

%Weight Vector
weight=w(theta_array);
%Log of Wieght Vector;

log_w=log(weight);
% histogram plot

 nbins = 25;
subplot(2,2,3);
hist(log_w,nbins);
subplot(2,2,4);
log_w(log_w<-0.3)=[];
hist(log_w,nbins);
% Expectation
Expect_2= sum(theta_array.*weight)/sum(weight)

% Varaince
Vari_2= sum((theta_array.^2).*weight)/sum(weight) - (Expect_2^2)


%Effective sample Size.
rev_weight=(weight.*S)/sum(weight);
S_effect=vpa(inv(sum(rev_weight.^2)))
