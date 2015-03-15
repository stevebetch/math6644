%Creates the A matrix for problem 5!
function [A,b,u]=SetupProb5b(n)
A=full(gallery('tridiag',n,1,-2,1));
h=1/n;
h2=h^2;
for i=1:n
    b(i,1)=h2*(((2*i)/n)-.5);
end
% b=zeros(n,1);

b(1)=1;
b(end)=-1;
u=0:1/(n-1):1; %creates a 1000x1 vector bounded by [0,1] 
end