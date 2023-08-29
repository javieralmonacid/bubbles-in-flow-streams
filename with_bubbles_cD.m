% With bubble problem, method of lines solution
%
% Author: Javier Almonacid
%         javiera@sfu.ca
%
% Date: June 15, 2023
%
clear; clc;

%w_in = 10;          % Inflow velocity [m/s]
%b = 1;              % Characteristic length 
%D = 1;              % Diffusion coefficient
%Pe = w_in * b / D;

%% Single Figure

Pe = 1e+03;
w = @(x) x.*(1-x);  % Non-dimensional velocity (may be zero along the boundary)
dx = 0.01;         % Spatial mesh size

alpha = 0.5;
beta = 50;
f = @(c) alpha*(1+tanh(beta*(c-1)));

[X,Z,C] = find_cD_not_shifted(dx,Pe,w,f);
figure(1)
surf(X,Z,C,'EdgeColor','none')
view([0 90])
max(max(C))
colorbar
colormap('jet')

xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
title('\bf Concentration $\widehat{c}_D$','Interpreter','latex','FontSize',14)

%% Contours of csat

figure(2)
clf
Pe = [250 500 750 1000 1250 1500 1750 2000];
w = @(x) x.*(1-x);  % Non-dimensional velocity (may be zero along the boundary)
dx = 0.01;         % Spatial mesh size
csat = 1e-03;

tic
for k = 1:length(Pe)
   [X,Z,C] = find_cD_not_shifted(dx,Pe(k),w,f);
   disp(['Solved for Pe =',num2str(Pe(k))])
   contour(X,Z,C,[csat; csat],'LineWidth',1.5)
   hold on
   grid on
end
toc

txt = '\leftarrow Pe = 250';
text(0.54,0.9,txt,'FontSize',12)

txt = '\leftarrow Pe = 500';
text(0.38,0.8,txt,'FontSize',12)

txt = '\leftarrow Pe = 750';
text(0.3,0.7,txt,'FontSize',12)

txt = '\leftarrow Pe = 1000';
text(0.25,0.6,txt,'FontSize',12)

txt = 'Pe = 2000 \rightarrow';
text(0.02,0.9,txt,'FontSize',12)

%%
xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
title('\bf Contours of $\widehat{c}_D = \widehat{c}_{sat} = \cdot 10^{-3}$','Interpreter','latex','FontSize',14)