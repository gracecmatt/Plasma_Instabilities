clear; clc; 

%%% ---------- change these values with parameters ----------
paramARR = linspace(0.25,0.75);
xaxis = "k"; 
filenameEnding = "k.png"; 
%%% ---------- end change section ----------
k = 0.5;
sigma = 1.5;
mu = 1;
syms omega D;
gamma = []; 
Omega = []; 
tic;

for i = 1:length(paramARR)
    k=paramARR(i); %%% ---------- change ----------

    % Dielectric function with kappa=2
    D(omega) = 1+(-1).*k.^(-2).*(2.*((-1).*k.*mu+omega).^2+k.^2.*sigma.^2).^(-3).*(8.* ...
          k.^2.*((-1).*k.*mu+omega).^4+24.*k.^4.*((-1).*k.*mu+omega).^2.*sigma.^2+ ...
          (sqrt(-1)*16).*2.^(1/2).*k.^5.*(k.*mu+(-1).*omega).*sigma.^3+(-6).* ...
          k.^6.*sigma.^4);

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
title("Kappa: All Roots \gamma("+xaxis+")",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #'); 
xlabel(xaxis); ylabel("\gamma("+xaxis+")");
grid on; colororder(newcolors)
% saveas(gcf,"Matlab/AllSolImages/gammaVs"+filenameEnding)

figure
plot(paramARR,Omega,'.-','linewidth',1)
title("Kappa: All Roots \Omega("+xaxis+")",'FontSize',14)
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