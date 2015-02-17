%Creates the A matrix for problem 5!

n=1000;
A=full(gallery('tridiag',n,1,-2-4*(1/1000),1));
b=zeros(n,1);
b(1)=1;
b(end)=-2;
u=0:1/(n-1):1; %creates a 1000x1 vector bounded by [0,1] 
