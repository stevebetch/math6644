function  [ x,count ] = GS( A,b,mindiff )
%GS Custom implementation of the Seidel method
%   This function takes in the A matrix, the b vector, the maximum number
%   of iterations, and the minimum error. returns the x values.

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
n=length(A);
while(1)
    count=count+1;
%     x=deinv*x0+dfinv;
%     for i=1:n
%         
%        x(i,1)=(1/D(i,i))*(b(i)-A(i,1:n)*x0+A(i,i)*x0(i));
%     end
    x=(b-A*x0+D*x0).*(1./d);
    if(max(abs(x-x0))<mindiff)
        break;
    end

    x0=x;
end


end



