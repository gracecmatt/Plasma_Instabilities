clear; clc

karr = 0.1:0.05:1.0;
theta1 = 1;
theta2 = 1;
mu1 = 0;
mu2 = 4;
beta = 0.9;
syms omega F;
gamma = []; 
Omega = []; 

for i = 1:length(karr)
    k = karr(i);
    F(omega) = 1+(sqrt(-1)*16).*2.^(1/2).*beta.*((1/4).*(2.^(1/2).*omega+(-1).*k.*(2.^( ...
        1/2).*mu1+(sqrt(-1)*3).*theta1)).*((sqrt(-1)*(-2)).*k.*mu1+(sqrt(-1)*2) ...
        .*omega+2.^(1/2).*k.*theta1).^(-3)+(-2).*k.^3.*(k.*mu1+(-1).*omega).* ...
        theta1.^3.*((-4).*k.*mu1.*omega+2.*omega.^2+k.^2.*(2.*mu1.^2+theta1.^2)) ...
        .^(-3))+(sqrt(-1)*(-16)).*2.^(1/2).*(1+(-1).*beta).*((1/4).*((sqrt(-1)*( ...
        -1)).*2.^(1/2).*(k.*mu2+(-1).*omega)+3.*k.*theta2).*(2.*k.*mu2+(-2).* ...
        omega+sqrt(-1).*2.^(1/2).*k.*theta2).^(-3)+(-2).*k.^3.*((-1).*k.*mu2+ ...
        omega).*theta2.^3.*(2.*((-1).*k.*mu2+omega).^2+k.^2.*theta2.^2).^(-3));
    
    Y = double(vpasolve(F(omega)==0,omega));
    gamma = [gamma, imag(Y)]; 
    Omega = [Omega, real(Y)]; 
end
numSol = length(gamma(:,1));

% Run to get only unique solutions
[gamma, Omega, numSol] = getUnique(gamma, Omega);

newcolors = {'#4363d8','#e6194B','#3cb44b','#42d4f4','#f032e6','#000075',...
    '#fabed4','#f58231','#469990','#9A6324','#800000','#aaffc3','#a9a9a9'};
colororder(newcolors)

plot(karr,gamma,'.-','linewidth',2)
title("Kappa Bump-on-Tail: All Roots \gamma(k)",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #');
xlabel("k"); ylabel("\gamma(k)");
grid on;

plot(karr,Omega,'.-','linewidth',2)
title("Kappa Bump-on-Tail: All Roots \Omega(k)",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #');
xlabel("k"); ylabel("\gamma(k)");
grid on;

function [gamma, Omega, numSol] = getUnique(gamma, Omega)
    klength = length(gamma(1,:));
    ImReMatrix = [gamma, Omega];

    % ----from MATLAB Answers user Walter Roberson 14 August 2020----
    [~,A] = uniquetol(ImReMatrix,'byrows',true); 
    ImReMatrix=ImReMatrix(sort(A),:); % unique rows in original order
    % ---------------------------------------------------------------

    gamma = ImReMatrix(:,1:klength);
    Omega = ImReMatrix(:,klength+1:end);
    numSol = length(gamma(:,1));
end