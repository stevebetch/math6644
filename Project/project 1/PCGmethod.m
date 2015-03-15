function [ x,count ] = PCGmethod( A,b,B,tol )
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
z=(r\B)'; %solve for z
p=z;

count=0;
debug=1;
 h=waitbar(0,'PCG method');
while debug>tol && count <100000
    count=count+1;
    if(mod(count,500)==0)
         waitbar(1,h,sprintf('%s%d','Running PCG! Iteration: ',count));
    end
    
    q=A*p;
    alpha=(z'*r)/(p'*q);
    
    x=x+alpha*p;
    ro=r;
    r=r-alpha*q;
    zo=z;
    z=(r\B)';%solve for z.
    beta=(z'*r)/(zo'*ro);
    
    p=z+beta*p;
    
    debug=(norm(r));
end
close(h);
end