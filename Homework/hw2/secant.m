function [ x, numIts  ] = secant( fhandle, x0,xneg,atol, rtol,maxIt)
%SECANT The secant method for non linear systems of equations
%%This method takes in the function handle to the system that needs to be
%%solved, the inital x value, and the tolerance.

r0=norm(fhandle(x0),inf);
x=x0;
x0=xneg;
fx=fhandle(x);
numIts=0;


while norm(fx,inf)>rtol*r0+atol && numIts<maxIt
    numIts=numIts+1;
    
  
    
    a=(fx-fhandle(x0))/(x-x0);
    
    x0=x;
    x=x0-fx/a;
    fx=fhandle(x);
    
    
end


end

