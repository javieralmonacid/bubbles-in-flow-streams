function out = find_concentration_dissolved_gas(reduced_pe, ...
                                                dx, ...
                                                Da_t, ...
                                                Da_s, ...
                                                lambda, ...
                                                xi)
% FIND_CONCENTRATION_DISSOLVED_GAS finds the concentration
% of dissolved gas in the system using the method of lines.
%
%  out = FIND_CONCENTRATION_DISSOLVED_GAS(reduced_pe, dx, Da_t) finds 
%  the concentration of dissolved gas in the system assuming no bubbles.
%  This assumes a constant flux of gas (Da_t) on the left boundary of 
%  the channel.
%
%  out = FIND_CONCENTRATION_DISSOLVED_GAS(reduced_pe, dx, Da_t, Da_s, 
%  lambda,xi) finds the concentration of dissolved gas in the system
%  with bubble formation.
%
%  Input parameters:
%    reduced_pe: reduced Peclet number
%    dx:         mesh spacing in x
%    Da_t:       Damkohler number for the generation of dissolved gas
%    Da_s:       Damkohler number for the separation of bubbles
%    lambda:     ratio between the rate of bubble generation and bubble
%                advection
%    xi:         dimensionless quantity that controls the rate of
%                bubble generation
%
%  Output parameters:
%    out (struct, nondimensional):
%      out.X: meshgrid for x coordinate
%      out.Z: meshgrid for z coordinate
%      out.C: concentration of dissolved gas.
%
%  Authors: Javier Almonacid, Simon Fraser University
%           Matthew Shirley, University of Oxford
%  Date:    June 15, 2023
%    

% Create mesh in x
x = 0:dx:1;
N = length(x);

if nargin == 3
    qleft = @(cd) Da_t + 0*cd;
elseif nargin == 6
    qs = @(cd) 1 + tanh(xi * (cd - 1));
    qleft = @(cd) (Da_t - Da_s * qs(cd))./(1 + Da_s * lambda * qs(cd));
else
    error('Incorrect number of arguments. Number of arguments allowed: 3 or 6')
end

% Mass matrix
M = diag([0, reduced_pe * x .* (1-x), 0]); % Mass matrix
options = odeset('Mass',M,'AbsTol',1e-8,'RelTol',1e-10,'MassSingular','yes');
c0 = zeros(N+2,1);
sol = ode15s(@(t,y) rhs_dissolved_gas_system(t,y,qleft,dx), ...
             [0,1], c0, options);

C = (deval(sol, x))';  % Evaluate function on uniform grid
C = C(:,2:end-1);      % Remove info from ghost nodes
[X,Z] = meshgrid(x,x); % Mesh for non-dimensional space variables.

out.X = X;
out.Z = Z;
out.C = C;


