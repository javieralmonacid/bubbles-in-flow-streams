% No bubbles problem, analytical solution

Pe = 100;
dx = 0.01;
dz = dx;

x = 0:dx:1;
z = 0:dz:1;
N = length(x); % No. discretization points in the x direction

[X,Z] = meshgrid(x,z);

% Solution in written in terms of series. We need to truncate it
% at some number 'cutoff'
cutoff = 200;
C = zeros(size(X));

for n = 1:cutoff
    C = C - 2/(n*pi)^2*( exp(-(n*pi)^2*Z/Pe) - 1).*cos(n*pi*X);
end

surf(X,Z,C,'EdgeColor','none')
colorbar
max(max(C))
view([0 90])