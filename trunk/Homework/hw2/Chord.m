function [x, numIts ] = Chord( fhandle, x0,atol, rtol,maxIt)
%CHORD The Chord method for non linear systems of equations
%%This method takes in the function handle to the system that needs to be
%%solved, the inital x value, and the tolerance.

r0=norm(fhandle(x0),inf);
x=x0;
fx=fhandle(x0);
numIts=0;
h=1e-5;

    for i=1:length(x)
        jacobian(i)=imag(fhandle(x+h*1i))/h;
    end
    [l,u]=lu(jacobian);
    

while norm(fx,inf)>rtol*r0+atol && numIts<maxIt
    numIts=numIts+1;
    
    if(mod(numIts,10000)==0)
        fprintf('%s %d\n','Chord Method: Max number of iterations:',numIts);
    end
   
    
    %% Solve for s. Since this is a 1-d problem, we can just devide by df.
    %if this was a multi-dimensional matrix, we would need to use a linear
    %solver to find s. 
    y=-fx/l;
    s=y/u;
    x=x+s;
    fx=fhandle(x);
    
    
end


end

