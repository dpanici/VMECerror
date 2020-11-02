% Run this script after running force_error.m to get the MHD energy
%% Energy
% W = W_magnetic + W_pressure
%W_magnetic = |B|^2 /2/mu0 , W_pressure = p / (gamma-1)
gamma=data.gamma;
presf = repmat(data.presf',1,dimU,dimV);
magB = sqrt((BU.^2).*dot(eu,eu,4) + (BV.^2).*dot(ev,ev,4));

[s1,u1,v1] = meshgrid(s,u,v);

% W_B = trapz(v,trapz(u,trapz(s,magB.*g))) / 2/mu0 
W_B = trapz(u,trapz(s,magB.*g)) / 2/mu0 *2*pi;



W_p = trapz(v,trapz(u,trapz(s,g.*presf./(gamma-1))));

W = W_B + W_p;
