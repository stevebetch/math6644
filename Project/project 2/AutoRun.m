%% This is the auto running script for Project 2

atol=1e-6;
rtol=1e-6;
maxIt=10000;
x0=ones(200,1);
fhandle=@(x) Eqn(x);
m=2;

Fstart=tic;
[ Fx, Fcount, Ferr ] = FixedPoint( fhandle, x0,atol, rtol,maxIt);
fend=toc(Fstart);

cstart=tic;
[Cx, Ccount, Cerr ] = Chord( fhandle, x0,atol, rtol,maxIt);
cend=toc(cstart);

nstart=tic;
[Nx,Ncount, Nerr]=Newton( fhandle, x0,atol, rtol,maxIt);
nend=toc(nstart);

shstart=tic;
[Shx,Shcount, Sherr]=Shamanskii( fhandle, x0,m,atol, rtol,maxIt);
shend=toc(shstart);





%% Ploting section figure 1
figure;

semilogy(0:Fcount-1,Ferr,'r--*','DisplayName','Fixed Point Error')
hold on; 
semilogy(0:Ccount-1,Cerr,'c--^','DisplayName','Chord Method Error');
semilogy(0:Ncount-1,Nerr,'b--d','DisplayName','Newton''s Method Error');
semilogy(0:Shcount-1,Sherr,'k--*','DisplayName','Shamanskii''s Method Error');


% Create ylabel
ylabel({'Relative error of the Chandrasekhar H-equation '},'FontSize',20,'Interpreter','latex');

% Create xlabel
xlabel({'Iteration Number'},'FontSize',20,'Interpreter','latex');

% Create title
title({'Error of Each Numerical Method at Each Iteration'},'FontSize',20,...
    'Interpreter','latex');

% Create legend
legend1 = legend('show');
set(legend1,'Interpreter','latex','FontSize',16);

hold off;



%% Latex out
fprintf('function  &Iterations Required to converge & Time to Converge & Time per Iteration &Final Error%s\n','\\\hline')
fprintf('Fixed Point&%d& %f&%f& %d%s \n',Fcount-1,fend,fend/(Fcount-1),Ferr(end),'\\\hline')
fprintf('Chord&%d& %f& %d&%d%s \n',Ccount-1,cend,cend/(Ccount-1),Cerr(end),'\\\hline')
fprintf('Newton''s Method&%d& %.8f&%d &%d%s\n',Ncount-1,nend,nend/(Ncount-1),Nerr(end),'\\\hline')
fprintf('Shamanskii''s Method&%d& %.8f&%d& %d%s \n',Shcount-1,shend,shend/(Shcount-1),Sherr(end),'\\\hline')
fprintf('\n');

