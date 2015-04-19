function [x, numIts ,stopCheck] = Chord( fhandle, x0,atol, rtol,maxIt)
%CHORD The Chord method for non linear systems of equations
%%This method takes in the function handle to the system that needs to be
%%solved, the inital x value, and the tolerance.

r0=norm(fhandle(x0),inf);
x=x0;
fx=fhandle(x0);
numIts=1;

    Jac=diffjac(x,fhandle,fx);
    [l,u]=lu(Jac);
    
stopVals=rtol*r0+atol;
stopCheck=r0;
while stopCheck(numIts)>stopVals && numIts<maxIt
   
    
    if(mod(numIts,10000)==0)
        fprintf('%s %d\n','Chord Method: Max number of iterations:',numIts);
    end
   
    
    %% Solve for s. Since this is a 1-d problem, we can just devide by df.
    %if this was a multi-dimensional matrix, we would need to use a linear
    %solver to find s. 
        [y,count1]=SOR(l,-fx,10^-6,1.1);
        [s,count2]=SOR(u,y,10^-6,.95);
        x=x+s;
    fx=fhandle(x);
     numIts=numIts+1;
    stopCheck(numIts)=norm(fx,inf);
    
end


end

