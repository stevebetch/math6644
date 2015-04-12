function [ Fx ] = Eqn( x )
%EQN This function takes in an x array and returns an Fx array. 
%   The c and N values are hard coded per the problem statement.
c=.9;
N=200;
u=(1:N)';
uj2=(u-.5)/N;

for i=1:length(x)
    ui=(i-.5)/N;
    insum=sum((ui*x)./(ui+uj2));
    Fx(i,1)=x(i)-(1-(c/(2*N)*insum))^-1;
end

end

