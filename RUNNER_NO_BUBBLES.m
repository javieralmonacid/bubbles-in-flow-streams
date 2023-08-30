% -----------------------------------------------------
% BUBBLES IN FLOW STREAMS AND POROUS MEDIA
% MPI 2023, June 12-16, 2023, 
% New Jersey Institute of Technology, Newark, NJ, USA
% -----------------------------------------------------
% 
% No bubbles problem, method of lines solution
%
% Author: Javier Almonacid
%         javiera@sfu.ca
%
% Date: August 29, 2023
%
clear; clc; close all;

%% Single computation

% Input parameters
reduced_Pe = 2000;  % Reduced Peclet number
dx = 0.01/4;        % Spatial mesh size
Da_t = 25;          % Damkohler number

% Compute non-dimensional solution
sol = find_concentration_dissolved_gas(reduced_Pe,dx,Da_t);

figure(1)
surf(sol.X,sol.Z,sol.C,'EdgeColor','none')
view([0 90])
max(max(sol.C))
colorbar
colormap('jet')

xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
title(['Concentration $\widehat{c}_d$ ($\mathcal{P} =$ ',num2str(reduced_Pe), ...
       ', $Da_t =$ ', num2str(Da_t), ')'], ...
       'Interpreter', 'latex', 'FontSize', 14)

%% Computation over 3 different Pe numbers

Pe = [500 1000 1500]; % Reduced Peclet numbers
dx = 0.01/4;          % Spatial mesh size
Da_t = 25;            % Damkohler number

max_cd = 0;

figure(2)
for i = 1:length(Pe)
    sol = find_concentration_dissolved_gas(Pe(i),dx,Da_t);
    if i == 1
        max_cd = max(max(sol.C));
    end
    
    subplot(1,length(Pe),i)
    surf(sol.X,sol.Z,sol.C,'EdgeColor','none')
    view([0 90])
    colorbar
    clim([0 max_cd])
    colormap('jet')
    xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
    ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
    title(['$\widehat{c}_d$ ($\mathcal{P} =$ ',num2str(Pe(i)), ...
           ', $Da_t =$ ', num2str(Da_t), ')'], ...
           'Interpreter', 'latex', 'FontSize', 14)
end
set(gcf,'Position',[977 726 1107 239])
exportgraphics(gcf,'figs/nobubble_01.pdf', 'ContentType', 'vector')

%% Computation over a series of (reduced) Pe numbers to compute cd=cd_sat

figure(3)
clf
Pe = [250 500 750 1000 1250 1500 1750 2000];
dx = 0.01; % Spatial mesh size
Da_t = 25;          % Damkohler number

tic
for k = 1:length(Pe)
   sol = find_concentration_dissolved_gas(Pe(k), dx, Da_t);
   % Rescale
   disp(['Solved for Reduced Pe = ',num2str(Pe(k))])
   contour(sol.X,sol.Z,sol.C,[1; 1],'LineWidth',1.5)
   hold on
   grid on
end
toc

txt = '$\leftarrow \mathcal{P} = 250$';
text(0.27,0.9,txt,'FontSize',14,'Interpreter','latex')

txt = '$\leftarrow \mathcal{P} = 500$';
text(0.18,0.8,txt,'FontSize',14,'Interpreter','latex')

txt = '$\leftarrow \mathcal{P} = 750$';
text(0.135,0.7,txt,'FontSize',14,'Interpreter','latex')

txt = '$\leftarrow \mathcal{P} = 1000$';
text(0.105,0.6,txt,'FontSize',14,'Interpreter','latex')

txt = '$\leftarrow \mathcal{P} = 2000$';
text(0.05,0.4,txt,'FontSize',14,'Interpreter','latex')

xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
title('\bf Contours of $c_d = c_{d,sat}$ - No bubbles','Interpreter','latex','FontSize',14)

exportgraphics(gcf,'figs/nobubble_02.pdf', 'ContentType', 'vector')