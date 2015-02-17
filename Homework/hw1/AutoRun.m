
SetupProb5;
h=waitbar(0,'1');
clear val;
w=1.49375;
for i=1:10
    waitbar(i/10,h,sprintf('%s%d%s','Running! ',i,' of 10 running!'));
val(i)=1/(10^i);
[xJ,countJ(i)]=Jacobi(A,b,val(i));
[xGS,countGS(i) ] = GS( A,b,val(i) );
[xSOR,countSOR2(i)]=SOR(A,b,val(i),w);

end
close(h);


figure1 = figure('Name','Convergence Criteria vs Number of iterations');

% Create axes
axes1 = axes('Parent',figure1,'YMinorTick','on','XMinorTick','on',...
    'XScale','log',...
    'XDir','reverse');
box(axes1,'on');
hold(axes1,'on');

% Create ylabel
ylabel('Number of iterations required to meet convergence criteria');

% Create xlabel
xlabel('Convergence Criteria');

% Create title
title('Convergence Criteria vs Number of iterations');

% Create semilogx
semilogx(val,countJ,'Color',[0 0 1]);
hold on;
semilogx(val,countGS,'g');
semilogx(val,countSOR2,'r');
set(semilogx1(1),'DisplayName','Gauss-Jacobi Method','Color',[0 0 1]);
set(semilogx1(2),'DisplayName','Gauss-Seidel Method','Color',[0 1 0]);
set(semilogx1(2),'DisplayName','SOR Method with \omega =1.49375','Color',[1 0 0]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',14);

clear countSOR;
control=20;
countSOR=zeros(10,1);
w=.2:(1.7-.2)/control:1.7;

h=waitbar(0,'1');
for i=1:control
waitbar(i/control,h,sprintf('%s%d%s%d%s','Running SOR ',i,' of ',control,' running!'));
    [xSOR,countSOR(i)]=SOR(A,b,10^-6,w(i));

end
close(h);
