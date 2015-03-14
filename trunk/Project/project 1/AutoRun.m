%THIS FUNCTION WILL CREATE MANY OF THE PLOTS SEEN IN THE PROJECT WRITE UP. RUN
%THIS SCRIPT TO SEE THE NUMERICAL METHODS IN ACTION.

clear;
ex=8:15;
runs=2.^ex+2;%+2 needed to make n=interior mesh points.
h=waitbar(0,'1');
clear val;
for a=1:length(runs)
    n=runs(a);
 waitbar(a/length(runs),h,sprintf('%s%d%s%d%s','Running! ',a,' of ', length(runs) ,' running!'));

[A,b,u]=SetupProb1(n);

S=zeros(runs(a));
Stic=tic;
k=1:runs(a);
parfor j=1:runs(a) 
        S(j,k)=sqrt(2/(n+1))*sin((pi*j*k)/(n+1));     
end
Ts(a)=toc(Stic(a));


PCGtic=tic;
A2=S*A*S';
[PCGx,PCGcount(a) ] = PCGmehtod( A2,b,(A2)^-1);
Tpcg(a)=toc(PCGtic(a));

CGtic=tic;
[CGx,CGcount(a)]=CGmehtod(A,b);
Tcg(a)=toc(CGtic(a));
end

close(h);


figure1 = figure('Name','Convergence Criteria vs Number of Iterations');

% Create axes
axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','YMinorTick','on',...
    'XMinorTick','on',...
    'XScale','log',...
    'XDir','reverse');
box(axes1,'on');
hold(axes1,'on');

% Create ylabel
ylabel('Number of iterations required to meet convergence criteria',...
    'FontSize',16,...
    'Interpreter','latex');

% Create xlabel
xlabel('Convergence Criteria','FontSize',16,'Interpreter','latex');

% Create title
title('Convergence Criteria vs Number of iterations','FontSize',16,'Interpreter','latex');

% Create semilogx
p1=semilogx(val,countJ,'Color',[0 0 1]);
set(p1,'DisplayName','Gauss-Jacobi Method','Color',[0 0 1]);
hold on;
p2=semilogx(val,countGS,'--g');
set(p2,'DisplayName','Gauss-Seidel Method','Color',[0 1 0]);
p3=semilogx(val,countSOR2,'--*r');


set(p3,'DisplayName','SOR Method with \omega =0.98','Color',[1 0 0]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',14);
hold off;
clear countSOR;
control=30;
countSOR=zeros(10,1);
w=.1:(1-.1)/(control-1):1;

h=waitbar(0,'1');
for i=1:control
waitbar(i/control,h,sprintf('%s%d%s%d%s','Running SOR ',i,' of ',control,' running!'));
    [xSOR,countSOR(i)]=SOR(A,b,10^-6,w(i));

end
figure2 = figure('Name','\omega values');
plot(w,countSOR);

close(h);
