%THIS FUNCTION WILL CREATE MANY OF THE PLOTS SEEN IN THE PROJECT WRITE UP. RUN
%THIS SCRIPT TO SEE THE NUMERICAL METHODS IN ACTION.

clear;
ex=8:14;
runs=2.^ex+2;%+2 needed to make n=interior mesh points.
h=waitbar(0,'1');
clear val;
for a=1:length(runs)
 waitbar(a/length(runs),h,sprintf('%s%d%s%d%s','Running! ',a,' of ', length(runs) ,' running!'));
 
 n=runs(a);
[A,b,u]=SetupProb1(n);

S=zeros(runs(a));
Stic(a)=tic;
k=1:n;
parfor j=1:n 
        S(j,k)=sqrt(2/(n+1))*sin((pi*j*k)/(n+1));     
end


Ts(a)=toc(Stic(a));
% 
% 
PCGtic(a)=tic;
A2=S*A*S';
[PCGx,PCGcount(a) ] = PCGmethod( A2,b,(A2)^-1,1e-6,0);
Tpcg(a)=toc(PCGtic(a));

CGtic(a)=tic;
[CGx,CGcount(a)]=CGmethod(A,b);
Tcg(a)=toc(CGtic(a));
end

close(h);

