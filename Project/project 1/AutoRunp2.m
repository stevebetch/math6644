%THIS FUNCTION WILL CREATE MANY OF THE PLOTS SEEN IN THE PROJECT WRITE UP. RUN
%THIS SCRIPT TO SEE THE NUMERICAL METHODS IN ACTION.

clear;
ex=0:6;
runs=50*2.^ex; %50,100,200,400,800...
p=[2,1,1/10,1/100];

h=waitbar(0,'1');
clear val;
for i=1:length(p)
    for a=1:length(runs)
     waitbar(a/length(runs),h,sprintf('%s%d%s%d%s%d','Running! ',a,' of ', length(runs) ,' running! P= ',p(i)));

     n=runs(a);
    [A,b,u]=SetupProb1(n,p(i));

    

%Begin with Strang's method
    TsS(a)=toc(Stic(a));
    PCGticS(a)=tic;
    [PCGxS,PCGcountS(a) ] = PCGmethod( A,b,C);
    TpcgS(a)=toc(PCGticS(a));

    CGtic(a)=tic;
    [CGxS,CGcount(a)]=CGmethod(A,b);
    TcgS(a)=toc(CGtic(a)); %#ok<*SAGROW>
    end
    
%  Now for Chan's method.

TsC(a)=toc(Stic(a));
    PCGticC(a)=tic;
    [PCGxC,PCGcountC(a) ] = PCGmethod( A,b,C);
    TpcgS(a)=toc(PCGticC(a));

    CGtic(a)=tic;
    [CGxC,CGcount(a)]=CGmethod(A,b);
    Tcg(a)=toc(CGtic(a));
    


end
close(h);

