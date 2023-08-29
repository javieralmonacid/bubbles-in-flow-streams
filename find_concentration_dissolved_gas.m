function out = find_concentration_dissolved_gas(reduced_pe, ...
                                                dx, ...
                                                Da_t, ...
                                                Da_s, ...
                                                lambda, ...
                                                xi)
% FIND_CONCENTRATION_DISSOLVED_GAS(reduced_pe, dx) finds the concentration
% of dissolved gas in the system using the method of lines.

% Create mesh in x
x = 0:dx:1;

if nargin == 2
    qleft = @(cd) 1 + 0*cd;
elseif nargin == 6
    qs = @(cd) 1 + tanh(xi * (cd - 1));
    qleft = @(cd) (Da_t - Da_s * qs(cd))./(1 + Da_s * lambda * qs(cd));
else
    error('Incorrect number of arguments. Number of arguments allowed: 2 or 6')
end

% Mass matrix
M = diag([0; reduced_pe * x .* (1-x); 0]); % Mass matrix
options = odeset('Mass',M,'AbsTol',1e-8,'RelTol',1e-10,'MassSingular','yes');
c0 = zeros(N+2,1);
sol = ode15s(@(t,y) rhs_dissolved_gas_system(t,y,qleft,dx), ...
             [0,1], c0, options);

% Evaluate function on uniform grid
[X,Z] = meshgrid(x,x);
C = (deval(sol, x))';


