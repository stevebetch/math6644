function [ x,count ] = CGmehtod( A,b,tol )
%CGMEHTOD Conjugate Gradient Method
%   This method takes in an A matrix and b vector as well as a final
%   tolerance. This method attempts to solve the system of linear Equations
%   A*x=b for x. The nxn matrix must be SPD and the b matrix should be of
%   length n. 



if( nargin==2) %only give the A/b matrix/vec
    tol=1e-8; %set the default tolerance to 10^-8
end

x=b; %Set x = b for the first guess.
r=b-A*x; %find the residual.

p=r;
if(norm(r)<tol)
    return;
end
count=0;
 h=waitbar(0,'CG method');
while (max(norm(r)))>tol
%     if(mod(count,100)==0)
%          waitbar(1,h,sprintf('%s%d','Running! Iteration: ',count));
%     end
    count=count+1;
    q=A*p;
  
    alpha=(r'*r)/(p'*q);
    x=x+alpha*p;
    ro=r;
    r=ro-alpha*q;
    beta=(r'*r)/(ro'*ro);
    p=r+beta*p;
   
end
close(h);