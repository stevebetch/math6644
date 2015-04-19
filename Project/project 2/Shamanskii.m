function [ x, numIts,stopCheck  ] = Shamanskii( fhandle, x0,m,atol, rtol,maxIt)
%Shamanskii Shamanskii iteration for non linear systems of equations
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
    [L,U]=lu(Jac);
    
    %% Solve for s. Since this is multi-dimensional matrix,
    %we need to use a linear solver to find s. 
    for a=1:m
        [y,count1]=SOR(L,-fx,10^-6,1.1);
        [s,count2]=SOR(U,y,10^-6,.95);
        x=x+s;
        fx=fhandle(x);
         numIts=numIts+1;
        stopCheck(numIts)=norm(fx,inf);
        if(stopCheck(numIts)<=stopVals)
            break;
        end
    end
end


end
