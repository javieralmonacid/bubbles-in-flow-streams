% No bubble problem, method of lines solution
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

Pe = 1e+02;

w = @(x) x.*(1-x);  % Non-dimensional velocity (may be zero along the boundary)
qleft = -1;         % Non-dimensional flux on left boundary.
dx = 0.005;          % Spatial mesh size

%% Main solver

x = 0:dx:1;
N = length(x); % No. discretization points in the x direction

% Mass matrix
M = diag([0; Pe*w(x'); 0]);
options = odeset('Mass',M,'AbsTol',1e-4,'RelTol',1e-6,'MassSingular','yes');

c0 = -x.*(1-x/2)*qleft;
c0 = [0; c0'; 0];
[z,C] = ode15s(@(t,y) dae_function(t,y,qleft,dx),[0,1],c0,options);
C = C(:,2:end-1); % Remove info from ghost nodes

[X,Z] = meshgrid(x,z);      % Mesh for non-dimensional space variables.
C = C + X.*(1-X./2).*qleft; % Non-dimensional concentration,
surf(X,Z,C,'EdgeColor','none')
view([0 90])
max(max(C))
colorbar

%xlabel('x')
%ylabel('z')
%title('Concentration of dissolved gas')

%% Contours of csat

csat = 5e-3;
contour(X,Z,C,[csat; csat])
set(gca,'XScale','log')




