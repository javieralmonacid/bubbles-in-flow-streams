function out = find_concentration_bubbles(reduced_fr, ...
                                          Da_s, ...
                                          lambda, ...
                                          xi, ...
                                          epsilon, ...
                                          sol_cd)
% FIND_CONCENTRATION_BUBBLES finds the concentration of bubbles
% in the system using the method of characteristics.
%
%  out = FIND_CONCENTRATION_BUBBLES(...) finds the concentration of
%  bubbles for a given concentration of dissolved gas.
%
%  Input parameters:
%    reduced_fr: reduced Froude number
%    dx:         mesh spacing in x
%    Da_s:       Damkohler number for the separation of bubbles
%    lambda:     ratio between the rate of bubble generation and bubble
%                advection
%    xi:         dimensionless quantity that controls the rate of
%                bubble generation
%    epsilon:    ratio of bubble radius to channel width
%    sol_cd:     structure coming from FIND_CONCENTRATION_DISSOLVED_GAS
%
%  Output parameters:
%    out (struct, nondimensional):
%      out.C: Concentration of bubbles that can be plot along sol_cd.X and
%      sol_cd.Z
%
%  Authors: Javier Almonacid, Simon Fraser University
%           Matthew Shirley, University of Oxford
%  Date:    June 15, 2023
%  

% Create function from data points along the left boundary of cd
z_left = sol_cd.Z(:,1);
c_left = sol_cd.C(:,1);
cd_z_fun = @(z) interp1(z_left,c_left,z);

% Define flux of bubbles
qs = @(cd) 1 + tanh(xi * (cd - 1));

% Solution can be computed using method of characteristics (see report)
gamma = @(z) Da_s * lambda * qs(cd_z_fun(z)) / ...
            (1 + Da_s * lambda * qs(cd_z_fun(z)));
CB = zeros(size(sol_cd.C));

% Fill entries of CB one by one
for i = 1:size(sol_cd.C,1)
    for j = 1:size(sol_cd.C,2)
        xx = sol_cd.X(i,j);
        zz = sol_cd.Z(i,j);
        aux = zz - (log(xx/(2*epsilon)) - 2*log((1-2*xx)/(1-4*epsilon)) + log((1-xx)/(1-2*epsilon))) / reduced_fr^2;
        if isreal(aux) && (aux > 0) && (aux < Inf)
            CB(i,j) = gamma(aux);
        else
            CB(i,j) = 0;
        end
    end
end

% Output concentration
out.C = CB;

% Output characteristic where CB = 0
xx = sol_cd.X(1,:);
zz = (log(xx/(2*epsilon)) - 2*log((1-2*xx)/(1-4*epsilon)) + log((1-xx)/(1-2*epsilon))) / reduced_fr^2;

out.char_x = xx;
out.char_z = real(zz);

