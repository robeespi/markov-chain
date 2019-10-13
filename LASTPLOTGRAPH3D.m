% Discrete Time Markov Chain
% Zombie Epidemic Model
% Transition Matrix,limiting distribution and Graph of Probability Distribution
%Cleaning the memory
clear all
%Set graphics object properties
set(gca,'FontSize',18);
set(0,'DefaultAxesFontSize',18);
%Setting the time steps
time=200;
%Seeting the size of the interval time for a DTMC
dtt=0.001; % Time step
%Instantiating infectious rate
beta=5*dtt;
%Instantiating recovery rate
gama=4*dtt;
% Instantiating population size
N=40; 
% plot every enth time interval
en=20; 
%Instantiating the transition matrix T
T=zeros(N+1,N+1); 
%Instantiating a vector for holding values of I
i=linspace(0,N,N+1);
%Instantiating a vector for holding probabilities for each time step 
p=zeros(time+1,N+1);
% Five individuals initially infected, the boat load.
p(1,2)=1; 
%Probability of increasing the number of infected
bt=(beta*i.*(N-i))/(N-1);
%Probability of decreasing in the number of infected
dt=i.*gama;
%Absorption state
T(1,1)=1;
%Fill a Element outside the loop
T(2,1)=dt(2);
%Loop to fill the transition matrix
for i=2:N % Define the transition matrix
T(i,i)=1-bt(i)-dt(i); % diagonal entries
T(i,i+1)=dt(i); % superdiagonal entries
T(i+1,i)=bt(i+1); % subdiagonal entries
end
%Filling elements last row of the transition matrix
T(N+1,N)=dt(N+1);
T(N+1,N+1)=1-dt(N+1);
%Running a loop to get a probability vector of a each time step 
for t=1:time
y=T*p(t,:)';
p(t+1,:)=y';
end
%Storing the initial distribution
pm(1,:)=p(1,:);
%Running a loop to store the probability of each time step
for t=1:time/en;
pm(t+1,:)=p(en*t,:);
end
%Instantiating a vector for holding the time steps for plotting
ti=linspace(0,time,time/en+1);
%Instantiating a vector for holding number of infectives for plotting
st=linspace(0,N,N+1);
%Draw a wireframe
mesh(st,ti,pm);
%Labeling 
xlabel('Number of Infectives');
ylabel('Time Steps');
zlabel('Probability');
%Setting the azimuth and elevation angles
view(140,30);
%Setting the axis limit
axis([0,N,0,time,0,1]);
%Code for get the stationary distribution of the transition matrix
A = T.'- eye(N+1);
A(N+1,:) = ones(1,N+1);
b = zeros(N+1,1);
b(N+1) = 1;
pi = A\b;
%Code for get cdf of I(t)
J=pi.*T;
TTT=ones(time+1,N+1);
YYY=ones-pm;





