function [ x, numIts  ] = FixedPoint( fhandle, x0,atol, rtol)
%FixedPoint Fixed-point iteration for non linear systems of equations
%%This method takes in the function handle to the system that needs to be
%%solved, the inital x value, and the tolerance.
%%Since the homework only requires a single nonlinear equation to be
%%sloved, this method was not extended to cover a matrix. 

r0=norm(fhandle(x0),inf);
x=x0;
fx=fhandle(x0);
numIts=0;
h=1e-5;

while norm(fx,inf)>rtol*r0+atol
    numIts=numIts+1;
    df=imag(fhandle(x+h*1i))/h;
    
    %% Remove this for final runtime calcs!
%     df2=(fhandle(x+h)-fx)/h;
%     dabs=abs(df-df2);
%     if(dabs>.001)
%         disp('ERROR WITH THE IMAGINARY STEP!!!');
%     end
    
    %% Solve for s. Since this is a 1-d problem, we can just devide by df.
    %if this was a multi-dimensional matrix, we would need to use a linear
    %solver to find s. 
    
    s=-fx/df;
    x=x+s;
    fx=fhandle(x);
    
    
end


end
