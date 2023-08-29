function out = rhs_dissolved_gas_system(~, y, f, dx)
% RHS_DISSOLVED_GAS_SYSTEM implements the right-hand side for the system
% that computed the concentration of dissolved gas in the systems

N = length(y)-2;
D2 = diag(-2*ones(N+2,1)) + diag(ones(N+1,1),1) + diag(ones(N+1,1),-1);
D2 = D2(2:end-1,:);

out = zeros(N+2,1);
out(1)       = (y(3) - y(1))/(2*dx) + f(y(2)); % Left flux-type BC
out(2:end-1) = (1/dx^2) * D2 * y;
out(end)     = (y(end)-y(end-2))/(2*dx);       % Right no-flux BC