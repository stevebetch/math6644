%% This is the auto running script for Project 2

atol=1e-6;
rtol=1e-6;
maxIt=10000;
%% Problem 2, part a

%Setup the problem 
f1= @(x) 2*x.^2-5; % create a function handle for the problem
x0=10; %initialize the x val 
xneg=0.99*x0; %for the secant method, x_(-1) val.
 fprintf('%s\n','Run 1');
[Nx1,newCount1]=Newton(f1,x0,atol,rtol);
[Cx1, Ccount1 ] = Chord( f1, x0,atol, rtol,maxIt);
[ Sx1, Scount1  ] = secant( f1, x0,xneg,atol, rtol,maxIt);
%% Problem 2, part b
 fprintf('%s\n','Run 2');

f2= @(x) sin(x)+x;
x0=0.5;
xneg=0.99*x0;

[Nx2,newCount2]=Newton(f2,x0,atol,rtol);
[Cx2, Ccount2 ] = Chord( f2, x0,atol, rtol,maxIt);
[ Sx2, Scount2  ] = secant( f2, x0,xneg,atol, rtol,maxIt);
%% Problem 2, part c

 fprintf('%s\n','Run 3');
f3= @(x) cos(x);
x0=03;
xneg=0.99*x0;


[Nx3,newCount3]=Newton(f3,x0,atol,rtol);
[Cx3, Ccount3 ] = Chord( f3, x0,atol, rtol,maxIt);
[ Sx3, Scount3  ] = secant( f3, x0,xneg,atol, rtol,maxIt);

%% Ploting section figure 1
figure;
x=-20:1:20;
plot(x,f1(x),'DisplayName','Function a')
hold on; 
plot(10,f1(10),'r*','DisplayName','Starting Point');
plot(Nx1,f1(Nx1),'c^','DisplayName','Newton''s method End Point');
plot(Cx1,f1(Nx1),'rd','DisplayName','Chord method End Point');
plot(Sx1,f1(Nx1),'k*','DisplayName','Secant method End Point');
plot(x,zeros(length(x),1),'g','DisplayName','Zero-line');

% Create ylabel
ylabel({'Y-Axis'},'FontSize',20,'Interpreter','latex');

% Create xlabel
xlabel({'X-Axis'},'FontSize',20,'Interpreter','latex');

% Create title
title({'Function a Start and End Points'},'FontSize',20,...
    'Interpreter','latex');

% Create legend
legend1 = legend('show');
set(legend1,'Interpreter','latex','FontSize',16);

hold off;
%% Ploting section Figure 2
figure;
x=-4:.1:4;
plot(x,f2(x),'DisplayName','Function b');

hold on;

plot(0.5,f2(0.5),'r*','DisplayName','Starting Point');
plot(Nx2,f2(Nx2),'c^','DisplayName','Newton''s method End Point');
plot(Cx2,f2(Cx2),'rd','DisplayName','Chord method End Point');
plot(Sx2,f2(Sx2),'k*','DisplayName','Secant method End Point');

plot(x,zeros(length(x),1),'g','DisplayName','Zero-line');
% Create ylabel
ylabel({'Y-Axis'},'FontSize',20,'Interpreter','latex');

% Create xlabel
xlabel({'X-Axis'},'FontSize',20,'Interpreter','latex');

% Create title
title({'Function b Start and End Points'},'FontSize',20,...
    'Interpreter','latex');

% Create legend
legend1 = legend('show');
set(legend1,'Interpreter','latex','FontSize',16);

hold off;
%% Ploting section Figure 3
figure;
x=-8:.1:4;
plot(x,f3(x),'DisplayName','Function c');

hold on;

plot(03,f3(3),'r*','DisplayName','Starting Point');
plot(Nx3,f3(Nx3),'r^','DisplayName','Newton''s method End Point');
% plot(Cx2,f3(Cx2),'rd');
% plot(Sx3,f3(Sx3),'k*');
plot(x,zeros(length(x),1),'g','DisplayName','Zero-line');
% Create ylabel
ylabel({'Y-Axis'},'FontSize',20,'Interpreter','latex');

% Create xlabel
xlabel({'X-Axis'},'FontSize',20,'Interpreter','latex');

% Create title
title({'Function c Start and End Points'},'FontSize',20,...
    'Interpreter','latex');

% Create legend
legend1 = legend('show');
set(legend1,'Interpreter','latex','FontSize',16);

hold off;


%% Latex out
fprintf('function  &Newton''s & Chord Method & Secant Method %s','\\\hline')
fprintf('a&%.8f& %.8f& %.8f%s \n',Nx1,Cx1,Sx1,'\\\hline')
fprintf('b&%.8f& %.8f& %.8f%s \n',Nx2,Cx2,Sx2,'\\\hline')
fprintf('c&%.8f& %s& %s%s \n',Nx3,'Did Not Converge','Did Not Converge','\\\hline')
fprintf('\n');

fprintf('function  &Newton''s & Chord Method & Secant Method %s\n','\\\hline')
fprintf('a&%d& %d& %d%s \n',newCount1,Ccount1,Scount1,'\\\hline')
fprintf('b&%d& %d& %d%s \n',newCount2,Ccount2,Scount2,'\\\hline')
fprintf('c&%d& %s& %s%s \n',newCount3,'Did Not Converge','Did Not Converge','\\\hline')
fprintf('\n');


