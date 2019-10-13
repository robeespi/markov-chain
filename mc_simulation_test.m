%  Markov chain simulation
clear all
mu=[0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];     % initial distribution
dtt=0.01; % Time step
beta=2*dtt;
gama=1*dtt;
k=100; % Total population size
P=zeros(k+1,k+1); % T is the transition matrix, defined below
i=linspace(0,k,k+1);
bt=(beta*i.*(k-i))/(k-1);
dt=i.*gama;
P(1,1)=1;
P(2,1)=dt(2);
for i=2:k % Define the transition matrix
P(i,i)=1-bt(i)-dt(i); % diagonal entries
P(i,i+1)=dt(i); % superdiagonal entries
P(i+1,i)=bt(i+1); % subdiagonal entries
end
P(k+1,k)=dt(k+1);
P(k+1,k+1)=1-dt(k+1);

n=20;           % number of steps to simulate
x=zeros(1,n+1); % instantiating x
t=0:n;          % time indexes

x(1)=rando(mu); % generate value of x(0) with distribucion mu

for i=1:n
  x(i+1) = rando(P(x(i),:));
end

g(1)=rando(mu);
for i=1:n
  g(i+1) = rando(P(x(i),:));
end

subplot(2,1,1)
stairs(t,x);
axis([0 n 0 (length(mu)+1)]);
title('A TRAJECTORY OF THE CHAIN');

subplot(2,1,2)
hist(x,length(mu))
title('Relative frequencies of the number of visits to states')
