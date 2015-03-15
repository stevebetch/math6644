%THIS FUNCTION WILL CREATE MANY OF THE PLOTS SEEN IN THE PROJECT WRITE UP. RUN
%THIS SCRIPT TO SEE THE NUMERICAL METHODS IN ACTION.

clear;
ex=0:8;
runs=50*2.^ex; %50,100,200,400,800...
p=[2,1,1/10,1/100];
% runs=4;
h=waitbar(0,'1');

for i=1:length(p)
    for a=1:length(runs)
     waitbar(a/length(runs),h,sprintf('%s%d%s%d%s%d','Running! ',a,' of ', length(runs) ,' running! P= ',p(i)));

     n=runs(a);
    [A,b]=SetupProb2a(n,p(i));

%%     Create Strang's Circulent matrix.
% CS=zeros(size(A));
CS=diag(diag(A)); % get the diagonal parts.
cfill=A(1,:); %grab first row. 
cfill(1)=[]; % kill off the diagonal part.
n2=floor(length(A)/2);
cfill=cfill(1:n2); % drop the excess data, 
cfill2=fliplr(cfill);
cfill2(1)=[];
cfill=[cfill,cfill2];
CS(1,2:end)=cfill;
        for j=2:length(A)
%             diagC=CS(j,j);
            cfillJ=circshift(cfill,j-1,2);
            CS(j,j+1:end)=cfillJ(j:end);
            CS(j,1:j-1)=cfillJ(1:j-1);
    
        end
%% Crete Chang's method

CC=diag(diag(A)); % get the diagonal parts.
Ccfill=A(1,:); %grab first row. 
% Ccfill(1)=[]; % kill off the diagonal part.
cj=zeros(1,length(A)-1);
for j=1:length(cj)
  cj(j)=((n-j)*Ccfill(j+1)+j*Ccfill(mod(-j,n+1)) )/n;
end

        parfor j=1:length(A)
            cjJ=circshift(cj,j-1,2);
            CC(j,j+1:end)=cjJ(j:end);
            CC(j,1:j-1)=cjJ(1:j-1);
        end


%% Create the DFT matrix.

% Ns=length(CS);
% Ws=ones(Ns);
% r=2:Ns;
% q=(2:Ns)';
% Ws(2:end,2:end)=exp(-(2.*pi.*q*r)/Ns)/sqrt(Ns);
% inWs=Ws^-1;




    
    
%Begin with Strang's method
    lambdaS=fft(CS(:,1));

    
    PCGticS(a,i)=tic;
%     a2A;
    [PCGxS,PCGcountS(a,i) ] = PCGmethod( A,b,CS^-1);
    TpcgS(a,i)=toc(PCGticS(a,i));

    CGtic(a,i)=tic;
    [CGxS,CGcount(a,i)]=CGmethod(A,b);
    TcgS(a,i)=toc(CGtic(a,i)); %#ok<*SAGROW>
    
    
%  Now for Chan's method.

    PCGticC(a,i)=tic;
%     a2=CC*A*CC';
    [PCGxC,PCGcountC(a,i) ] = PCGmethod( A,b,CC^-1);
    TpcgS(a,i)=toc(PCGticC(a,i));

    CGtic(a,i)=tic;
    [CGxC,CGcount(a,i)]=CGmethod(A,b);
    Tcg(a,i)=toc(CGtic(a,i));
    

    end
end
close(h);
%% NEXT SECTION: PART B
%% 
%% 
%%
%%

h=waitbar(0,'1');
for i=1:length(p)
    for a=1:length(runs)
     waitbar(a/length(runs),h,sprintf('%s%d%s%d%s%d','Running! ',a,' of ', length(runs) ,' running! P= ',p(i)));

     n=runs(a);
    [A,b]=SetupProb2b(n,p(i));
%     n=2*n;
%%     Create Strang's Circulent matrix.
% CSb=zeros(size(A));
CSb=diag(diag(A)); % get the diagonal parts.
Bcfill=A(1,:); %grab first row. 
Bcfill(1)=[]; % kill off the diagonal part.
Bn2=floor(length(A)/2);
Bcfill=Bcfill(1:Bn2); % drop the excess data, 
Bcfill2=fliplr(Bcfill);
Bcfill2(1)=[];
Bcfill=[Bcfill,Bcfill2];
CSb(1,2:end)=Bcfill;
        parfor j=2:length(A)
%             diagC=CSb(j,j);
            BcfillJ=circshift(Bcfill,j-1,2);
            CSb(j,j+1:end)=BcfillJ(j:end);
            CSb(j,1:j-1)=BcfillJ(1:j-1);
    
        end
%% Crete Chang's method

CCb=diag(diag(A)); % get the diagonal parts.
Ccfillb=A(1,:); %grab first row. 
% Ccfill(1)=[]; % kill off the diagonal part.
cjb=zeros(1,length(A)-1);
for j=1:length(cjb)
  cjb(j)=((n-j)*Ccfillb(j+1)+j*Ccfillb(mod(-j,n+1)) )/n;
end

        parfor j=1:length(A)
            cjJb=circshift(cjb,j-1,2);
            CCb(j,j+1:end)=cjJb(j:end);
            CCb(j,1:j-1)=cjJb(1:j-1);
        end


%% Create the DFT matrix.

% Ns=length(CSb);
% Ws=ones(Ns);
% r=2:Ns;
% q=(2:Ns)';
% Ws(2:end,2:end)=exp(-(2.*pi.*q*r)/Ns)/sqrt(Ns);
% inWs=Ws^-1;




    
    
%Begin with Strang's method
    lambdaS=fft(CSb(:,1));

    
    PCGticSB(a,i)=tic;
%     a2A;
    [PCGxSB,PCGcountSB(a,i) ] = PCGmethod( A,b,CSb^-1);
    TpcgSB(a,i)=toc(PCGticSB(a,i));

    CGticB(a,i)=tic;
    [CGxSB,CGcountB(a,i)]=CGmethod(A,b);
    TcgSB(a,i)=toc(CGticB(a,i)); %#ok<*SAGROW>
    
    
%  Now for Chan's method.

    PCGticCB(a,i)=tic;
%     a2=CC*A*CC';
    [PCGxCB,PCGcountCB(a,i) ] = PCGmethod( A,b,CCb^-1);
    TpcgSB(a,i)=toc(PCGticCB(a,i));

    CGBtic(a,i)=tic;
    [CGxCB,CGcount(a,i)]=CGmethod(A,b);
    TcgB(a,i)=toc(CGticB(a,i));
    

    end
end
close(h);


