function [X,Z,C,Sol] = find_cD_not_shifted(dx,Pe,w,f)
    x = 0:dx:1;
    N = length(x);

    M = diag([0; Pe*w(x'); 0]); % Mass matrix
    options = odeset('Mass',M,'AbsTol',1e-8,'RelTol',1e-10,'MassSingular','yes');
    c0 = zeros(N+2,1);
    
    Sol = ode15s(@(t,y) dae_function(t,y,f,dx),[0,1],c0,options);
    %z = (Sol.x)'; %Use ode15s points
    z = linspace(0,1,150); %equal spaced points
    %C = (Sol.y)'; %use ode15s
    C = (deval(Sol,z))'; %Equally spaced points
    
    %[z,C] = ode15s(@(t,y) dae_function(t,y,f,dx),[0,1],c0,options);

    [X,Z] = meshgrid(x,z);      % Mesh for non-dimensional space variables.
    C = C(:,2:end-1);           % Remove info from ghost nodes
end

function out = dae_function(~,y,f,dx)
    N = length(y)-2;
    D2 = diag(-2*ones(N+2,1)) + diag(ones(N+1,1),1) + diag(ones(N+1,1),-1);
    D2 = D2(2:end-1,:);

    out = zeros(N+2,1);
    out(1) = y(3)-y(1)+2*dx*(1-f(y(2)));
    out(2:end-1) = (1/dx^2)*D2*y;
    out(end) = (y(end)-y(end-2))/(2*dx);
end
