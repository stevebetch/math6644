%Creates the A matrix for problem 5!
function [A,b,u]=SetupProb1(n)
A=full(gallery('tridiag',n,1,-2-(1/n)^2*(2),1));
b=zeros(n,1);
b(1)=1;
b(end)=-1;
u=0:1/(n-1):1; %creates a 1000x1 vector bounded by [0,1] 
end