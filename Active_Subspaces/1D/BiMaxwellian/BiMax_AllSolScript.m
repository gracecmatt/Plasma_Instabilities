clear; clc

%%% ---------- change these values with parameters ----------
paramARR = 0.1:0.01:2; 
xaxis = "k"; 
filenameEnding = "k.png"; 
%%% ---------- end change ----------
k = 1;
sigma1 = 1;
sigma2 = 1;
mu1 = 0;
mu2 = 4;
beta1 = 0.9;  beta2 = 1-beta1;
syms z F;
gamma = []; 
Omega = []; 
tic;


% A1 = @(x) 1/sigma1*(x - mu1);
% A2 = @(x) 1/sigma2*(x - mu2);

% F(z) 1+(2*beta1/(sigma1^2*k^2))*(1+A1(z)*zetaf(A1(z)))+(2*beta2/(sigma2^2*k^2))*(1+A2(z)*zetaf(A2(z))) ;

% options = optimset('Display','off');
% z_r1 = fsolve(F, [guess_real, guess_imag], options);
% xi = abs(z_r1(1)) + 1i*z_r1(2); % complex-valued angular frequency ==(omega/k)



for i = 1:length(paramARR)
    k = paramARR(i); %%% ---------- change ----------

    % A1 = @(x) 1/sigma1*(x - mu1);
    % A2 = @(x) 1/sigma2*(x - mu2);

    A1 = @(zr, zi) 1/sigma1*(zr + 1i*zi- mu1);
    A2 = @(zr, zi) 1/sigma2*(zr + 1i*zi- mu2);

    % F(omega) = 1+(-2).*beta.*(2.*((-1).*k.*mu1+omega).^2+k.^2.*sigma1.^2).^(-3).*(4.*(( ...
    %       -1).*k.*mu1+omega).^4+12.*k.^2.*((-1).*k.*mu1+omega).^2.*sigma1.^2+(-3) ...
    %       .*k.^4.*sigma1.^4+(sqrt(-1)*8).*2.^(1/2).*k.^4.*mu1.*(sigma1.^2).^(3/2)+ ...
    %       (sqrt(-1)*(-8)).*2.^(1/2).*k.^3.*omega.*(sigma1.^2).^(3/2))+(-2).*((-1)+ ...
    %       beta).*(2.*((-1).*k.*mu2+omega).^2+k.^2.*sigma2.^2).^(-3).*((-4).*((-1) ...
    %       .*k.*mu2+omega).^4+(-12).*k.^2.*((-1).*k.*mu2+omega).^2.*sigma2.^2+3.* ...
    %       k.^4.*sigma2.^4+(sqrt(-1)*(-8)).*2.^(1/2).*k.^4.*mu2.*(sigma2.^2).^(3/2) ...
    %       +(sqrt(-1)*8).*2.^(1/2).*k.^3.*omega.*(sigma2.^2).^(3/2));
    
    % F(z) = 1+(2*beta1/(sigma1^2*k^2))*(1+A1(z)*zetaf(A1(z)))+(2*beta2/(sigma2^2*k^2))*(1+A2(z)*zetaf(A2(z))) ;
    F = @(z) 1+(2*beta1/(sigma1^2*k^2))*(1+A1(real(z),imag(z))*zetaf(A1(real(z),imag(z))))+(2*beta2/(sigma2^2*k^2))*(1+A2(real(z),imag(z))*zetaf(A2(real(z),imag(z))));

    Y = double(vpasolve(F(z)==0,z));

    gamma = [gamma, imag(Y)];
    Omega = [Omega, real(Y)];
end
numSol = length(gamma(:,1));

% Run to get only unique solutions
[gamma, Omega, numSol] = getUnique(gamma, Omega);

% Run to get only actual solutions
[gamma, Omega, numSol] = getActualSol(gamma, Omega, F, paramARR);


newcolors = {'#4363d8','#e6194B','#3cb44b','#42d4f4','#f032e6','#000075',...
    '#fabed4','#f58231','#469990','#9A6324','#800000','#aaffc3','#a9a9a9'};
colororder(newcolors)

figure(1)
plot(paramARR,gamma,'.-','linewidth',2)
title("Kappa Bump-on-Tail: All Roots \gamma("+xaxis+")",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #');
xlabel(xaxis); ylabel("\gamma("+xaxis+")");
grid on;
saveas(gcf,"Matlab/AllSolImages/gammaVs"+filenameEnding)

figure(2)
plot(paramARR,Omega,'.-','linewidth',2)
title("Kappa Bump-on-Tail: All Roots \Omega("+xaxis+")",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #');
xlabel(xaxis); ylabel("\Omega("+xaxis+")");
grid on;
colororder(newcolors)
saveas(gcf,"Matlab/AllSolImages/OmegaVs"+filenameEnding)
toc

function [gamma, Omega, numSol] = getActualSol(gamma, Omega, F, paramARR)
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
        Fzmax = max( double(abs(F(gamma(j,:)+1i*Omega(j,:)))) );
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