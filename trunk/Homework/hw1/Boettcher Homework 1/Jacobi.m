function [ x,count ] = Jacobi( A,b,mindiff )
%JACOBI Custom implementation of the jacobi method
%   This function takes in the A matrix, the b vector, the maximum number
%   of iterations, and the minimum error. returns the x values.
if(~any(diag(A)))
    error 'There is a diagonal entry that is 0!'
end

count=0;
x0=zeros(size(b));
while(1)
    count=count+1;
    x=(b-A*x0+(diag(A).*x0))./diag(A);
    
    if(abs(x-x0)<mindiff)
        break;
    end

    x0=x;
end


end


