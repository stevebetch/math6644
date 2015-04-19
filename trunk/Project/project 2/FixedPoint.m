function [ x, numIts,stopCheck ] = FixedPoint( fhandle, x0,atol, rtol,maxIt)
%FixedPoint Fixed-point iteration for non linear systems of equations
%%This method takes in the function handle to the system that needs to be
%%solved, the inital x value, and the tolerance.
%%Since the homework only requires a single nonlinear equation to be
%%sloved, this method was not extended to cover a matrix. 

r0=norm(fhandle(x0),inf);
x=x0;
fx=fhandle(x0);
numIts=1;


stopVals=rtol*r0+atol;
stopCheck=r0;

while stopCheck(numIts)>stopVals && numIts<maxIt
    
    x=x-fx;
    
    fx=fhandle(x);
    
        numIts=numIts+1;
    stopCheck(numIts)=norm(fx,inf);

    
end


end
