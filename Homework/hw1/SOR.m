function  [ x,count ] = SOR( A,b,mindiff,w )
%GS Custom implementation of the Seidel method
%   This function takes in the A matrix, the b vector, the minimum error 
%   and the relaxation parameter, w. returns the x values.

if(~any(diag(A)))
    error 'There is a diagonal entry that is 0!'
end

D=diag(diag(A));
d=diag(A);
E=triu(A-D);
F=tril(A-D);
deinv=((D-E)^-1)*F;
dfinv=(D-F)^-1*b;

count=0;
x0=zeros(size(b));
while(1)
    count=count+1;
    xgs=(b-A*x0+D*x0).*(1./d);
    x=w*xgs+(1-w)*x0;

    if((abs(x-x0)<mindiff))

        break;
    end
    
    if(count>100000)
       break;
    end

    x0=x;
end


end