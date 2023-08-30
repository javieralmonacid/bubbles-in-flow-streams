function [X,Z,C] = find_cD(dx,Pe,w)

x = 0:dx:1;
qleft = -1;

% Mass matrix
M = diag([0; Pe*w(x'); 0]);
options = odeset('Mass',M,'AbsTol',1e-8,'RelTol',1e-10,'MassSingular','yes');

c0 = -x.*(1-x/2)*qleft;
c0 = [0; c0'; 0];
[z,C] = ode15s(@(t,y) dae_function(t,y,qleft,dx),[0,1],c0,options);
C = C(:,2:end-1); % Remove info from ghost nodes

[X,Z] = meshgrid(x,z);      % Mesh for non-dimensional space variables.
C = C + X.*(1-X./2).*qleft; % Non-dimensional concentration (undo shifting)

end

function out = dae_function(t,y,qleft,dx)

N = length(y)-2;
D2 = diag(-2*ones(N+2,1)) + diag(ones(N+1,1),1) + diag(ones(N+1,1),-1);
D2 = D2(2:end-1,:);

out = zeros(N+2,1);
out(1) = (y(3)-y(1))/(2*dx);
out(2:end-1) = (1/dx^2)*D2*y - qleft;
out(end) = (y(end)-y(end-2))/(2*dx);

end
