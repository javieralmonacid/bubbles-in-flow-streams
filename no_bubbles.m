% No bubble problem, method of lines solution
clear all;

Pe = 10;
dx = 0.01;
qleft = -1;

x = 0:dx:1;
N = length(x); % No. discretization points in the x direction

M = diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
% Neumann correction
M(1,2)   = 2;
M(N,N-1) = 2;
M = 1/(Pe*dx^2)*M;

odefun = @(t,y) M*y - (1/Pe)*qleft*ones(length(y),1);
c0 = -x.*(1-x/2)*qleft;
c0 = c0';
[z,C] = ode23s(odefun,[0,1],c0);


[X,Z] = meshgrid(x,z);
C = C + X.*(1-X./2).*qleft;
surf(X,Z,C,'EdgeColor','none')
view([0 90])
max(max(C))
colorbar


