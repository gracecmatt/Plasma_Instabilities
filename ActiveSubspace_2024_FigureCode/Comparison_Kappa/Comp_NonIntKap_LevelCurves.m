%% Kappa Distribution Comparison Level Curves
% GM 1/11/24

clear; clc; 
%%% --------- change these values with parameter for level curve ---------
% paramARR = linspace(0.001,1,200);
paramARR = logspace(-1,0,102);
xaxis = "k"; 
filenameEnding = "k.png"; 
%%% ---------- end change section ----------
kappaIntARR = [2,3];
% kapARR = [2,2.2,2.45,3];
kapARR = [1.63,1.9,2.45,3.14];
% initialize vectors
omegaXieShiftScaled = zeros(length(kapARR),length(paramARR));
omegaExact = zeros(length(kapARR),length(paramARR));

% choose parameters
k = 0.5;
sigma = 1;
mu = 1;

% initialize for symbolic root-finder vpasolve()
syms omega D;

% %% find the "actual" solutions -> polynomial roots
% for j = 1:length(kappaIntARR)
%     kappa = kappaIntARR(j);
%     gamma = []; Omega = [];
%     for i = 1:length(paramARR)
%         k=paramARR(i); %%% -- level curve parameter ----------
% 
%         % dielectric functions exported from Wolfram Mathematica (w/ zero mean)
%         if (kappa==2)
%             D(omega) = 1+(-2).*k.^(-2).*(2.*omega.^2+k.^2.*sigma.^2).^(-3).*(4.*k.^2.*omega.^4+ ...
%               12.*k.^4.*omega.^2.*sigma.^2+(sqrt(-1)*(-8)).*2.^(1/2).*k.^5.*omega.* ...
%               sigma.^3+(-3).*k.^6.*sigma.^4);
%         elseif (kappa==3)
%             D(omega) = 1+(-2).*k.^(-2).*(2.*omega.^2+3.*k.^2.*sigma.^2).^(-4).*(8.*k.^2.* ...
%               omega.^6+60.*k.^4.*omega.^4.*sigma.^2+270.*k.^6.*omega.^2.*sigma.^4+( ...
%               sqrt(-1)*(-144)).*6.^(1/2).*k.^7.*omega.*sigma.^5+(-135).*k.^8.* ...
%               sigma.^6); 
%         elseif (kappa==4)
%             D(omega) = 1+(-2).*k.^(-2).*(2.*omega.^2+5.*k.^2.*sigma.^2).^(-5).*(16.*k.^2.* ...
%               omega.^8+224.*k.^4.*omega.^6.*sigma.^2+1400.*k.^6.*omega.^4.*sigma.^4+ ...
%               7000.*k.^8.*omega.^2.*sigma.^6+(sqrt(-1)*(-3200)).*10.^(1/2).*k.^9.* ...
%               omega.*sigma.^7+(-4375).*k.^10.*sigma.^8);
%         else
%             % error statement here
%         end
% 
%         % symbolic root finder finds all roots, Y is a complex vector
%         % mu*k is added, applying invariance of dielectric function
%         Y = double(vpasolve(D(omega)==0,omega)) + mu*k;
% 
%         gamma = [gamma, imag(Y)];
%         Omega = [Omega, real(Y)];
%     end
%     numSol = length(gamma(:,1));
% 
%     % Run to get only unique solutions
%     [gamma, Omega, numSol] = getUnique(gamma, Omega);
% 
%     % Run to get only actual solutions (real part is not zero or increasing for large k)
%     [gamma, Omega, numSol] = getActualSol(gamma, Omega, D, paramARR);
% 
%     for i = 1:length(paramARR)
%         k = paramARR(i);
% 
%         % find the root with the largest imaginary part
%         Y = Omega(:,i)+1i*gamma(:,i);
%         temp = sort(1i*Y,'ComparisonMethod','real');
%         omegaExact(j,i) = -1i*temp(1);
%     end
% end

%% compute solutions using Xie/Weideman algorithm

for i = 1:length(paramARR)
    k = paramARR(i);
    
    for j = 1:length(kapARR)
        kappa = kapARR(j);

        % Xie/Weideman algorithm
        init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k,sigma,0,kappa) + mu*k;
        xi_guess_shiftscaled = (init_guess/k-mu)/sigma; % shifted & scaled
        omegaXieShiftScaled(j,i) = Kappa_Disp_Using_Xie(k*sigma,1,0,kappa,xi_guess_shiftscaled)*sigma*k + mu*k;
    
        omegaExact(j,i) = Kappa_exact(k*sigma,kappa,xi_guess_shiftscaled)*sigma*k + mu*k;

        % omegaError(j,i) = abs(real(omegaXieShiftScaled(j,i))-real(omegaExact(j,i)))+1i*(imag(omegaXieShiftScaled(j,i))-imag(omegaExact(j,i)));
        omegaRelError(j,i) = abs(omegaExact(j,i)-omegaXieShiftScaled(j,i))/abs(omegaExact(j,i));
    end
   
end


%% plot figures 
% % with ideal formatting (as of 1/27/24)
close all;

% color = {'#dc143c','#cb2888','#9932cc','#6c4ce2','#2f5ada'};%red->blue
% color = {'#000000','#dc143c','#cb2888','#9932cc','#000000','#6c4ce2','#2f5ada'};
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure
% loglog(paramARR,imag(omegaExact(1,:)),'-','linewidth',2); hold on
% loglog(paramARR,imag(omegaXieShiftScaled(1,:)),'s','linewidth',2,'MarkerSize',5); 
% loglog(paramARR,imag(omegaXieShiftScaled(2,:)),'d','linewidth',2,'MarkerSize',4); 
% loglog(paramARR,imag(omegaXieShiftScaled(3,:)),'x','linewidth',2,'MarkerSize',6); 
% loglog(paramARR,imag(omegaExact(2,:)),'-','linewidth',2); 
% loglog(paramARR,imag(omegaXieShiftScaled(end,:)),'+','linewidth',2,'MarkerSize',6); 
% title('\textbf{Kappa $\gamma(k)$ Level Curves}','FontSize',16);
%     legend('Exact $(\kappa\!=\!2)$','$\kappa\!=\!2$','$\kappa\!=\!2.2$','$\kappa\!=\!2.45$','Exact $(\kappa\!=\!3)$','$\kappa\!=\!3$','Location','NorthEast','FontSize',12);
loglog(paramARR, omegaRelError,'linewidth',2);
title('Relative Error in $\gamma(k)$','FontSize',16)
    legend('$\kappa=1.63$','$\kappa=1.9$','$\kappa=2.45$','$\kappa=3.14$') %1.63,1.9,2.45,3.14
    xlabel('$k$','FontSize',14);
    % ylabel('$\gamma(k)$','FontSize',14);
    ylabel('$\frac{\left|\gamma-\gamma_{i} \right|}{|\gamma|}$','FontSize',14')
    set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
    % set(gca,'Yscale','log','XScale','Linear')
    % colororder(color)
    set(gcf, 'PaperPosition', [0 0 6 4.5]); %Position the plot further to the upper-left conder
    set(gcf, 'PaperSize', [6 4.5]); % Extends the plot to fill the entire paper
%     saveas(gcf, 'Figs/gammak_LevelCurves_nonIntKap_2.pdf')
%     saveas(gcf, 'Figs/gammak_LevelCurves_nonIntKap_2.fig')
%     hold off
% 
%     save('Data/Comparison_NonInt_kappa_data.mat')



%% functions
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