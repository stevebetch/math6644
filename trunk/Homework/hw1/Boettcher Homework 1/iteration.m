function [ out ] = iteration( A )
%ITERATION Summary of this function goes here
%   Detailed explanation goes here

itlimit=100;
a=0;

M=diag(diag(A));
N=M-A;

G=M^-1*N;
x0=[1;0;01];
c0=[0;1;0];
b=[3;4;5];
while (1)
    a=a+1;
    x1=G*x0+c0;
    
    diff=abs(sum(x1-x0)/3);
    
    if(diff<0.000001 || a>itlimit)
        
       break;
    end
    x0=x1;
    c0=M^-1*b;
end
out=x1;
end

