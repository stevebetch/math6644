function [ x,count ] = PCGmehtod( A,b,B,tol )
%CGMEHTOD Preconditioned Conjugate Gradient Method
%   This method takes in an A matrix and b vector as well as a final
%   tolerance. This method attempts to solve the system of linear Equations
%   A*x=b for x. The nxn matrix must be SPD and the b matrix should be of
%   length n. 

if( nargin==3) %only give the A/b matrix/vec
    tol=1e-8; %set the default tolerance to 10^-8
end

x=b; %Set x = b for the first guess.
r=b-A*x; %find the residual.
z=(r\B)';
p=z;

count=0;
debug=1;
while debug>tol && count <100000
    count=count+1;
    q=A*p;
    alpha=(z'*r)/(p'*q);
    
    x=x+alpha*p;
    ro=r;
    r=r-alpha*q;
    zo=z;
    z=(r\B)';
    beta=(z'*r)/(zo'*ro);
    
    p=z+beta*p;
    
    debug=(norm(r));
end