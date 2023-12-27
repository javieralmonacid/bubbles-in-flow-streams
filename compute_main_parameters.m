% -----------------------------------------------------
% BUBBLES IN FLOW STREAMS AND POROUS MEDIA
% MPI 2023, June 12-16, 2023, 
% New Jersey Institute of Technology, Newark, NJ, USA
% -----------------------------------------------------
% 
% Compute main parameters for the bubble problem
% (Pcal, Fcal, Da_t, Da_s, lambda, xi, epsilon)
%
% Author: Javier Almonacid
%         javiera@sfu.ca
%
% Date: August 30, 2023
%
function out = compute_main_parameters(delta_p, ...
                                       qt, ...
                                       Rcal_s, ...
                                       bubble_radius)

% Compute [w], hereafter called wmax
%delta_p = 2500;     % Pa (guessed)
h = 1;               % m
a = 1e-4;            % m
mu0_water = 8.9e-04; % Pa s
wmax = delta_p * a^2 / (8*h*mu0_water);

% Compute Pcal
delta = a/h;
D_hydrogen = 5.11e-09; % m^2/s
Pcal = delta^2 * wmax * h / D_hydrogen;

% Compute Fcal
g = 9.81; % m/s^2
Fcal = wmax / (delta * sqrt(g*h));

% Compute Da_t
%qt = 1e-02; % Guessed
cd_sat_hydrogen = 7.9e-4; % mol/L
Da_t = delta * h * qt / (D_hydrogen * cd_sat_hydrogen);

% Compute Da_s
%Rcal_s = 1e-10; % Guessed
Da_s = delta * h * Rcal_s / (D_hydrogen * cd_sat_hydrogen);

% Compute epsilon;
%bubble_radius = 1e-6; % Guessed
epsilon = bubble_radius / a;

% Compute lambda
Pe = Pcal/delta^2;
rho0_water = 997.05;      % kg/m^3
Re = rho0_water * wmax * h / mu0_water;
cb_max = cd_sat_hydrogen; % Guessed
mub_hydrogen = 9.03e-06;  % kg/(m s)
lambda = (3/2) * (1/(delta^2*Re)) * (cd_sat_hydrogen/cb_max) * ...
         (1/(2*epsilon^3*(1-2*epsilon)*(1-4*epsilon))) * ...
         ((2+3*mub_hydrogen/mu0_water)/(1-mub_hydrogen/mu0_water));

% Compute xi
beta = 1/cd_sat_hydrogen; % Guessed
xi = beta*cd_sat_hydrogen;

% Compute Fr
Fr = delta * Fcal;

% Display results
table(wmax,Pcal,Fcal,Da_t,Da_s,lambda,xi,epsilon,Re,Pe,Fr)

if nargout > 0
    out.reduced_Pe = Pcal;
    out.reduced_Fr = Fcal;
    out.Da_t = Da_t;
    out.Da_s = Da_s;
    out.lambda = lambda;
    out.xi = xi;
    out.epsilon = epsilon; 
end