clear; clc



%% Kappa
% bulk energy noralized to 1: E0 = 1/2*me*mu0^2 = 1
% we normalize all velocities thermal velocity
clear; clc
me = 9.109383701528*10^(-31); % [kg]

kBT = 2*10^5 * 8.617333262*10^(-5);
v_th = sqrt(kBT/me); % [m/s]

mu0 = sqrt(2/me) / v_th
sigma = (sqrt(2)*v_th) / v_th  % [1]

%% Incomplete Max
clear; clc;
me = 9.109383701528*10^(-31); % [kg]
kBT = 60; % [eV]
v_th1 = sqrt(kBT/me); % [m/s]

mu = sqrt(2/me) / v_th1; % [m/s]
nu = mu*0.5-1; % [m/s] ****
% mu = mu0/mu0-1; nu = nu0/mu0-1; v_th1 = v_th1/mu0-1; % [1] normalize by the mean

v_th2 = @(sig) sqrt(sig.^2/2 + sig.*(mu-nu).*exp(-(mu-nu)^2./sig.^2)./(sqrt(pi*sig.^2)/2.*(1+erf((mu-nu)./sig))) );
sigma_invf = @(sig) abs(v_th2(sig) - v_th1);

sigma = fsolve(sigma_invf,1)