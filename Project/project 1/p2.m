ex=0:6;
runs=50*2.^ex; %50,100,200,400,800...
p=1;
% runs=4;
h=waitbar(0,'1');
i=1;
    for a=1:length(runs)
     waitbar(a/length(runs),h,sprintf('%s%d%s%d%s','Running! ',a,' of ', length(runs) ,' running!'));

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
diagsB=diag(A);
        parfor j=2:length(A)
%             diagC=CSb(j,j);
            BcfillJ=circshift(Bcfill,j-1,2);
            CSb(j,:)=[BcfillJ(1:j-1),diagsB(j),BcfillJ(j:end)];
        end
%% Crete Chang's method

CCb=diag(diag(A)); % get the diagonal parts.
Ccfillb=A(1,:); %grab first row. 
% Ccfill(1)=[]; % kill off the diagonal part.
cjb=zeros(1,length(A)-1);
parfor j=1:length(cjb)
  cjb(j)=((n-j)*Ccfillb(j+1)+j*Ccfillb(mod(-j,n+1)) )/n;
end
diagb=diag(A);
        parfor j=1:length(A)
            cjJb=circshift(cjb,j-1,2);
            CCb(j,:)=[cjJb(1:j-1),diagb(j),cjJb(j:end)];
        end


        
%Begin with Strang's method
      
     
     a1=ifft(fft(CSb(:,1)).^-1)';
     a1s=zeros(length(a1));
     parfor j=1:length(A)
            a1s(:,j)=circshift(a1,j-1,2);
     end
     a0=a1s*A*a1s';
     a2=(a0)^-1;
     PCGticSB(a,i)=tic;
    [PCGxSB,PCGcountSB(a,i) ] = PCGmethod( a0,b,a2,1e-6,0);
    TpcgSB(a,i)=toc(PCGticSB(a,i));

    CGticB(a,i)=tic;
    [CGxSB,CGcountB(a,i)]=CGmethod(A,b);
    TcgSB(a,i)=toc(CGticB(a,i)); %#ok<*SAGROW>
    
    
%  Now for Chan's method.

    
     a1=ifft(fft(CCb(:,1)).^-1)';
     a1s=zeros(length(a1));
     parfor j=1:length(A)
            a1s(:,j)=circshift(a1,j-1,2);
     end
     a0=a1s*A*a1s';
     a2=(a0)^-1;
     PCGticCB(a,i)=tic;
    [PCGxCB,PCGcountCB(a,i) ] = PCGmethod( a0,b,a2,1e-6,0);
    TpcgCB(a,i)=toc(PCGticCB(a,i));

%     CGBtic(a,i)=tic;
%     [CGxCB,CGcount(a,i)]=CGmethod(A,b);
%     TcgCB(a,i)=toc(CGticB(a,i));
%     

    end

close(h);

