clear; clc; 

%%% ---------- change these values with parameters ----------
paramARR = linspace(0.15,0.55,250);
xaxis = "k"; 
filenameEnding = "k.pdf"; 
%%% ---------- end change section ----------
omegaAApprox = zeros(1,length(paramARR));
omegaXieShiftScaled = zeros(1,length(paramARR));
errorXie = zeros(1,length(paramARR));
errorAApprox = zeros(1,length(paramARR));
% k = 0.3;
kappa = 2;
sigma1 = 0.55;
sigma2 = 0.45;
mu = 10;
v0 = 4;
beta = 0.95;
syms omega D;
gamma = []; 
Omega = []; 
tic;

for i = 1:length(paramARR)
    k = paramARR(i); %%% ---------- change ----------

    % Dielectric function with kappa=2 (using mu=0, need to add mu*k later)
    D = @(omega) 1+(-1).*k.^(-2).*(2.*beta.*k.^2.*(2.*omega.^2+k.^2.*sigma1.^2).^(-3).*( ...
          4.*omega.^4+12.*k.^2.*omega.^2.*sigma1.^2+(sqrt(-1)*(-8)).*2.^(1/2).* ...
          k.^3.*omega.*sigma1.^3+(-3).*k.^4.*sigma1.^4)+(1/2).*(2.*((-1)+beta).* ...
          k.^2.*(k.^2.*sigma2.^2+2.*(omega+(-1).*k.*v0).^2).^(-3).*(3.*k.^4.* ...
          sigma2.^4+(sqrt(-1)*8).*2.^(1/2).*k.^3.*sigma2.^3.*(omega+(-1).*k.*v0)+( ...
          -12).*k.^2.*sigma2.^2.*(omega+(-1).*k.*v0).^2+(-4).*(omega+(-1).*k.*v0) ...
          .^4)+2.*((-1)+beta).*k.^2.*(k.^2.*sigma2.^2+2.*(omega+k.*v0).^2).^(-3).* ...
          (3.*k.^4.*sigma2.^4+(sqrt(-1)*8).*2.^(1/2).*k.^3.*sigma2.^3.*(omega+k.* ...
          v0)+(-12).*k.^2.*sigma2.^2.*(omega+k.*v0).^2+(-4).*(omega+k.*v0).^4))); 

    Y = double(vpasolve(D(omega)==0,omega));

    gamma = [gamma, imag(Y)];
    Omega = [Omega, real(Y) + mu*k]; % adding mu back in here
end
numSol = length(gamma(:,1));

% Run to get only unique solutions
[gamma, Omega, numSol] = getUnique(gamma, Omega);

% Run to get only actual solutions
[gamma, Omega, numSol] = getActualSol(gamma, Omega, D, paramARR);

%% Compute the analytic approximation and Xie/Weideman root
for i = 1:length(paramARR)
    k = paramARR(i);

    Y = Omega(:,i)+1i*gamma(:,i);
    temp = sort(1i*Y,'ComparisonMethod','real');
    exact = -1i*temp(1);

    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k,sigma1,sigma2,0,v0,beta,kappa) + mu*k;
    xi_guess_shiftscaled = (init_guess/k-mu)/sigma1;
    omegaXieShiftScaled(i) = BiKappa_Disp_Using_Xie(k*sigma1,1,sigma2/sigma1,0,v0/sigma1,beta,kappa,xi_guess_shiftscaled)*sigma1*k + mu*k;
    omegaAApprox(i) = AnalyticApproxKappa(k,sigma1,sigma2,v0,beta,kappa) + mu*k;
    
    errorXie(i) = abs(imag(exact)-imag(omegaXieShiftScaled(i)))./abs(imag(exact));
    errorAApprox(i) = abs(imag(exact)-imag(omegaAApprox(i)))./abs(imag(exact));
end
toc;
newcolors = {'#4363d8','#e6194B','#3cb44b','#42d4f4','#000075',...
    '#fabed4','#f58231','#469990','#9A6324','#800000','#aaffc3','#f032e6','#a9a9a9'};

%% Run Figures
close all;
% imaginary part of roots
figure
plot(paramARR,gamma,'.-','linewidth',1); hold on
plot(paramARR,imag(omegaXieShiftScaled),'dg','linewidth',2);
plot(paramARR,imag(omegaAApprox),'om','linewidth',1);
title("Kappa Bump-on-Tail: All Roots \gamma("+xaxis+")",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #'); 
xlabel(xaxis); ylabel("\gamma("+xaxis+")");
grid on; colororder(newcolors)
% saveas(gcf,"Matlab/AllSolImages/gammaVs"+filenameEnding)

% % real part of roots
% figure
% plot(paramARR,Omega,'.-','linewidth',1); hold on
% plot(paramARR,real(xie_shiftscaled),'dg','linewidth',2);
% plot(paramARR,real(omegaAApprox),'om','linewidth',1);
% title("Kappa Bump-on-Tail: All Roots \Omega("+xaxis+")",'FontSize',14)
% lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
% title(lgd,'Solution #');
% xlabel(xaxis); ylabel("\Omega("+xaxis+")");
% grid on; colororder(newcolors)
% % saveas(gcf,"Matlab/AllSolImages/OmegaVs"+filenameEnding)

% relative error between maximal root and Weideman/Xie and Analytic Approx
figure
semilogy(paramARR,errorXie,'.-','linewidth',1); hold on
semilogy(paramARR,errorAApprox,'.-','linewidth',1); 
title("Kappa Bump-on-Tail Relative Error in \gamma("+xaxis+")",'FontSize',14)
legend('Weideman/Xie Algorithm','Analytic Approximation','Location','Best');
xlabel(xaxis); ylabel("\gamma("+xaxis+") (log scale)");
grid on; colororder(newcolors)
set(gcf, 'PaperPosition', [0 0 6 5]); %Position the plot further to the upper-left conder
set(gcf, 'PaperSize', [6 5]); % Extends the plot to fill the entire paper
saveas(gcf, "Figs/gammaRelError_"+filenameEnding, 'pdf')

%% Functions
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