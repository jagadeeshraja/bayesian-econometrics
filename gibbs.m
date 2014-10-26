% Given parameters are as:
y=[1 0.5];
rho=0.8;
f1=zeros(11002,1);
f2=zeros(11002,1);
% Intial conditions are as:
f2(1)=0;
f1(1)=0;

% Sampling algorithm
for i=1:11000
   f1(i+1)=normrnd(y(1)+rho*(f2(i)-y(2)),1-rho^2);
   f2(i+1)=normrnd(y(2)+rho*(f1(i+1)-y(1)),1-rho^2);
end
% the required value of (x,y) with burn in 500 
x=f1(100:11000);
y=f2(100:11000);

% scatter plot
scatter(x,y);
