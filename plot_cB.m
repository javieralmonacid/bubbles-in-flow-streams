clear; clc; close all;

Pe = 100;
dx = 0.01;
dz = dx;
w = @(x) x.*(1-x);  % Non-dimensional velocity (may be zero along the boundary)

alpha = 0.5;
Beta = 75;


a = 1e-2;
b = 1;
delta = a/b;
R = 1e-07;
mu = 1e-03;
rho0 = 1e+03;
ndw = 1e-01;
ndcB = 1;
ndqL = 1e-10;
ndc = 1;
csat = 1e-4;

Rl = (3/2) * (delta/(R^3*(1-R)*(1-2*R))) ...
           * mu/(rho0 * b * ndw)...
           * ndqL/(ndw * ndcB);
       
Rl = 1;
       
% Get cD
qL_fun = @(cc) 1+tanh(Beta*ndc*(cc-csat/ndc)); 
[X,Z,CD,Sol] = find_cD_not_shifted(dx,Pe,w,qL_fun);

figure()
surf(X,Z,CD,'EdgeColor','none')
colorbar
view([0 90])
xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
title(['\bf Concentration $\widehat{c}_D$',', Pe = ',num2str(Pe)],'Interpreter','latex','FontSize',14)

aux = @(t) 0.05*(log(t/R) - 2*log((1-2*t)/(1-2*R)) + log((1-t)/(1-R)));

%cd_bdy = CD(:,1);
%z_bdy = Z(:,1);
%plot(z_bdy,cd_bdy)


for i=1:size(X,1)
    for j=1:size(X,2)
        if X(i,j) < 1/2
            if Z(i,j) < aux(X(i,j))
                CD(i,j) = 0;
            end
        else
            CD(i,j) = 0;
        end
    end
end


qL = qL_fun(CD);
CB = Rl * qL;




figure()
surf(X,Z,CB,'EdgeColor','none')
colorbar
view([0 90])
xlabel('$\hat{x}$','Interpreter','latex','FontSize',16)
ylabel('$\hat{z}$','Interpreter','latex','FontSize',16)
title(['\bf Concentration $\widehat{c}_B$',', Pe = ',num2str(Pe)],'Interpreter','latex','FontSize',14)





       
