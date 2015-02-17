function  [ x,count ] = GS( A,b,mindiff )
%GS Custom implementation of the Seidel method
%   This function takes in the A matrix, the b vector, the maximum number
%   of iterations, and the minimum error. returns the x values.

if(~any(diag(A)))
    error 'There is a diagonal entry that is 0!'
end

D=diag(diag(A));
E=triu(A-D);
F=tril(A-D);
deinv=((D-E)^-1)*F;
dfinv=(D-F)^-1*b;
count=0;
x0=zeros(size(b));
while(1)
    count=count+1;
    x=deinv*x0+dfinv;
    
    if(abs(x-x0)<mindiff)
        break;
    end

    x0=x;
end


end



