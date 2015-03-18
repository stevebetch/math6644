%THIS FUNCTION WILL CREATE MANY OF THE PLOTS SEEN IN THE PROJECT WRITE UP. RUN
%THIS SCRIPT TO SEE THE NUMERICAL METHODS IN ACTION.

clear;
ex=0:7;
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
diags=diag(A);
        parfor j=2:length(A)
%             diagC=CS(j,j);
            cjJ=circshift(cfill,j-1,2);
            CS(j,:)=[cjJ(1:j-1),diags(j),cjJ(j:end)];
        end
%% Crete Chang's method

CC=diag(diag(A)); % get the diagonal parts.
Ccfill=A(1,:); %grab first row. 
% Ccfill(1)=[]; % kill off the diagonal part.
cj=zeros(1,length(A)-1);
parfor j=1:length(cj)
  cj(j)=((n-j)*Ccfill(j+1)+j*Ccfill(mod(-j,n+1)) )/n;
end
diags=diag(A);
        parfor j=1:length(A)
            cjJ=circshift(cj,j-1,2);
            CC(j,:)=[cjJ(1:j-1),diags(j),cjJ(j:end)];
            
        end






    
    
%Begin with Strang's method


     a1=ifft(fft(CS(:,1)).^-1)';
     a1s=zeros(length(a1));
     parfor j=1:length(A)
            a1s(:,j)=circshift(a1,j-1,2);
     end

    PCGticS(a,i)=tic;
    [PCGxS,PCGcountS(a,i) ] = PCGmethod( A,b,a1s);
    TpcgS(a,i)=toc(PCGticS(a,i));

%     CGtic(a,i)=tic;
%     [CGxS,CGcount(a,i)]=CGmethod(A,b);
%     TcgS(a,i)=toc(CGtic(a,i)); %#ok<*SAGROW>
    
    
%  Now for Chan's method.

    
     a1=ifft(fft(CC(:,1)).^-1)';
     a1s=zeros(length(a1));
     parfor j=1:length(A)
            a1s(:,j)=circshift(a1,j-1,2);
     end
    PCGticC(a,i)=tic;
    [PCGxC,PCGcountC(a,i) ] = PCGmethod( A,b,a1s);
    TpcgC(a,i)=toc(PCGticC(a,i));

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
p2(); %done for debugging purposes. 
