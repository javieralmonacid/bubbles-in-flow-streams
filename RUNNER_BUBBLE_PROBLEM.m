% -----------------------------------------------------
% BUBBLES IN FLOW STREAMS AND POROUS MEDIA
% MPI 2023, June 12-16, 2023, 
% New Jersey Institute of Technology, Newark, NJ, USA
% -----------------------------------------------------
% 
% Bubble problem. Solution computed using method of 
% lines + method of characteristics.
%
% Author: Javier Almonacid
%         javiera@sfu.ca
%
% Date: August 29, 2023
%
clear; clc; close all;

%% Single computation

% Input parameters
reduced_Pe = 500;  % Reduced Peclet number
reduced_Fr = 4;
dx = 0.01/4;        % Spatial mesh size
Da_t = 500;          % Damkohler number (dissolved gas)
Da_s = 1000;          % Damkohler number (bubbles)
lambda = 1;
xi = 15;
epsilon = 1e-3;

% Compute non-dimensional solution
sol_CD = find_concentration_dissolved_gas(reduced_Pe,dx,Da_t,Da_s,lambda,xi);
max_CD = max(max(sol_CD.C));
min_CD = min(min(sol_CD.C));

% Compute concentration of bubbles
sol_CB = find_concentration_bubbles(reduced_Fr, Da_s, lambda, xi, epsilon, sol_CD);
max_CB = max(max(sol_CB.C));
min_CB = min(min(sol_CB.C));

%% Plot nondimensional solutions

min_colorbar = min([min_CD, min_CB]);
max_colorbar = max([max_CD, max_CB]);

figure(1)

subplot(1,2,1)
surf(sol_CD.X,sol_CD.Z,sol_CD.C,'EdgeColor','none')
view([0 90])
colorbar
clim([min_colorbar,max_colorbar])
colormap('jet')

xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
title(['$\widehat{c}_d$ ($\mathcal{P} =$ ',num2str(reduced_Pe), ...
       ', $Da_t =$ ', num2str(Da_t), ')'], ...
       'Interpreter', 'latex', 'FontSize', 14)

subplot(1,2,2)
surf(sol_CD.X,sol_CD.Z,sol_CB.C,'EdgeColor','none')
view([0 90])
colorbar
clim([min_colorbar,max_colorbar])
colormap('jet')
xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
title(['$\widehat{c}_b$ ($\mathcal{F} =$ ',num2str(reduced_Fr), ...
       ', $Da_s =$ ', num2str(Da_s), ')'], ...
       'Interpreter', 'latex', 'FontSize', 14)


