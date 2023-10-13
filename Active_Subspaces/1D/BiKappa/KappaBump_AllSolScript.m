clear; clc; 

%%% ---------- change these values with parameters ----------
paramARR = 0.2:0.01:1.5;
xaxis = "theta1"; 
filenameEnding = "theta1.png"; 
%%% ---------- end change section ----------
k = 0.9;
theta1 = 0.5;
theta2 = 1;
mu1 = 0;
mu2 = 4;
beta = 0.5;
syms omega D;
gamma = []; 
Omega = []; 
tic;

for i = 1:length(paramARR)
    theta1 = paramARR(i); %%% ---------- change ----------

    % D(omega) = 1+(-2).*beta.*(2.*((-1).*k.*mu1+omega).^2+k.^2.*theta1.^2).^(-3).*(4.*(( ...
    %       -1).*k.*mu1+omega).^4+12.*k.^2.*((-1).*k.*mu1+omega).^2.*theta1.^2+(-3) ...
    %       .*k.^4.*theta1.^4+(sqrt(-1)*8).*2.^(1/2).*k.^4.*mu1.*(theta1.^2).^(3/2)+ ...
    %       (sqrt(-1)*(-8)).*2.^(1/2).*k.^3.*omega.*(theta1.^2).^(3/2))+(-2).*((-1)+ ...
    %       beta).*(2.*((-1).*k.*mu2+omega).^2+k.^2.*theta2.^2).^(-3).*((-4).*((-1) ...
    %       .*k.*mu2+omega).^4+(-12).*k.^2.*((-1).*k.*mu2+omega).^2.*theta2.^2+3.* ...
    %       k.^4.*theta2.^4+(sqrt(-1)*(-8)).*2.^(1/2).*k.^4.*mu2.*(theta2.^2).^(3/2) ...
    %       +(sqrt(-1)*8).*2.^(1/2).*k.^3.*omega.*(theta2.^2).^(3/2));

    % Dielectric function with kappa=2
    D(omega) = 1+(sqrt(-1)*16).*2.^(1/2).*beta.*((1/4).*(2.^(1/2).*omega+(-1).*k.*(2.^( ...
      1/2).*mu1+(sqrt(-1)*3).*theta1)).*((sqrt(-1)*(-2)).*k.*mu1+(sqrt(-1)*2) ...
      .*omega+2.^(1/2).*k.*theta1).^(-3)+(-2).*k.^3.*(k.*mu1+(-1).*omega).* ...
      theta1.^3.*((-4).*k.*mu1.*omega+2.*omega.^2+k.^2.*(2.*mu1.^2+theta1.^2)) ...
      .^(-3))+(sqrt(-1)*16).*2.^(1/2).*(1+(-1).*beta).*((1/4).*(2.^(1/2).* ...
      omega+(-1).*k.*(2.^(1/2).*mu2+(sqrt(-1)*3).*theta2)).*((sqrt(-1)*(-2)).* ...
      k.*mu2+(sqrt(-1)*2).*omega+2.^(1/2).*k.*theta2).^(-3)+(-2).*k.^3.*(k.* ...
      mu2+(-1).*omega).*theta2.^3.*((-4).*k.*mu2.*omega+2.*omega.^2+k.^2.*(2.* ...
      mu2.^2+theta2.^2)).^(-3));

    Y = double(vpasolve(D(omega)==0,omega));

    gamma = [gamma, imag(Y)];
    Omega = [Omega, real(Y)];
end
numSol = length(gamma(:,1));

% Run to get only unique solutions
[gamma, Omega, numSol] = getUnique(gamma, Omega);

% Run to get only actual solutions
[gamma, Omega, numSol] = getActualSol(gamma, Omega, D, paramARR);


newcolors = {'#4363d8','#e6194B','#3cb44b','#42d4f4','#f032e6','#000075',...
    '#fabed4','#f58231','#469990','#9A6324','#800000','#aaffc3','#a9a9a9'};

 
figure
plot(paramARR,gamma,'.-','linewidth',1)
title("Kappa Bump-on-Tail: All Roots \gamma("+xaxis+")",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #'); 
xlabel(xaxis); ylabel("\gamma("+xaxis+")");
grid on; colororder(newcolors)
% saveas(gcf,"Matlab/AllSolImages/gammaVs"+filenameEnding)

figure
plot(paramARR,Omega,'.-','linewidth',1)
title("Kappa Bump-on-Tail: All Roots \Omega("+xaxis+")",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #');
xlabel(xaxis); ylabel("\Omega("+xaxis+")");
grid on; colororder(newcolors)
% saveas(gcf,"Matlab/AllSolImages/OmegaVs"+filenameEnding)
toc

function [gamma, Omega, numSol] = getActualSol(gamma, Omega, D, paramARR)
% this function does 2 checks:
% it evaluates |F(gamma+i*Omega)| < tol
% and makes sure the derivative of gamma(param) is small and not zero 

    solArr = [];
    tol = 5; % tolerance on solution value (about 0)
    dtol = 0.25; % tolerance on derivative (should be small or negative)
    colMax = length(gamma(1,:));
    numSol = length(gamma(:,1));
    z90 = floor(colMax*0.90); % index at 90% of length_k
    z88 = floor(colMax*0.88); % index at 88% of length_k
    for j=1:numSol
        Fzmax = max( double(abs(D(gamma(j,:)+1i*Omega(j,:)))) );
        dgammaEnd = (gamma(j,z90)-gamma(j,z88))/(paramARR(z90)-paramARR(z88)); %slope of gamma
        if Fzmax < tol && dgammaEnd < dtol && dgammaEnd ~= 0.0
            solArr = [solArr; j];
        end
    end
    
    gamma = gamma(solArr,:);
    Omega = Omega(solArr,:);
    numSol = length(gamma(:,1));
end


function [gamma, Omega, numSol] = getUnique(gamma, Omega)
    N = 3; % round to account for machine error giving extra duplicates
    ImReMatrix = [round(gamma,N), round(Omega,N)];

    % ----from MATLAB Answers user Walter Roberson 14 August 2020----
    [~,A] = uniquetol(ImReMatrix,'byrows',true); 
    ImReMatrix(sort(A),:); % unique rows in original order
    % ---------------------------------------------------------------

    gamma = gamma(sort(A),:);
    Omega = Omega(sort(A),:);
    numSol = length(gamma(:,1));
end