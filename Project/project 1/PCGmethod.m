function [ x,count ] = PCGmethod( A,b,B,tol )
%CGMEHTOD Preconditioned Conjugate Gradient Method
%   This method takes in an A matrix and b vector as well as a final
%   tolerance. This method attempts to solve the system of linear Equations
%   A*x=b for x. The nxn matrix must be SPD and the b matrix should be of
%   length n. 

if( nargin==3) %only give the A/b matrix/vec
    tol=1e-6; %set the default tolerance to 10^-8
end

x=b; %Set x = b for the first guess.
r=b-A*x; %find the residual.
% z=CGmethod(B,r);
z=(r\B)'; %solve for z
p=z;
r0=norm(r);
count=0;
debugV=r0;


while debugV>tol && count <100000
    count=count+1;
    if(mod(count,500)==0)
        if(~exist('h','var'))
         h=waitbar(0,'PCG method');
        end
         waitbar(1,h,sprintf('%s%d','Running PCG! Iteration: ',count));
    end
    
    q=A*p;
    alpha=(z'*r)/(p'*q);
    
    x=x+alpha*p;
    ro=r;
    r=r-alpha*q;
    zo=z;
%     z=CGmethod(B,r);
    z=(r\B)';%solve for z.
    beta=(z'*r)/(zo'*ro);
    
    p=z+beta*p;
    
   debugV=norm(r)/r0;
end
if(exist('h','var'))
    
    close(h);
end
end