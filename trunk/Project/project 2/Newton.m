function [ x, numIts ,stopCheck ] = Newton( fhandle, x0,atol, rtol,maxIt)
%NEWTON Newton's method for non linear systems of equations
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
   
    Jac=diffjac(x,fhandle,fx);
    
%% Solve for s. Since this is multi-dimensional matrix,
    %we need to use a linear solver to find s. 
    
    [s,count2]=SOR(Jac,-fx,10^-6,.95);
    x=x+s;
    fx=fhandle(x);
     numIts=numIts+1;
    stopCheck(numIts)=norm(fx,inf);
    
end


end
