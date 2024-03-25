%% 1D Kappa Bump-on-Tail All Solution Script Computing Relative Error in gamma
% between analytic approximation in Sarkar et al. (2015) and our numeric
% algorithm
clear; clc;  
%%% --- parameters to change when tuning the numerical algorithm ---
% L (zetaph) = 3
% M (spectral) = 2^14
% Vmax (spectral) = 200

%%% ---------- change these values to plot other level curves ----------
% paramARR = logspace(-0.78,0,100);
paramARR = logspace(-1,0,100);
kapARR = [2,4];

xaxis = "k"; 
filenameEnding = "k.pdf"; 
%%% ---------- end change section ----------
omegaExact = zeros(length(paramARR),length(kapARR));
omegaAApprox = zeros(length(paramARR),length(kapARR));
omegaXieShiftScaled = zeros(length(paramARR),length(kapARR));
errorXie = zeros(length(paramARR),length(kapARR));
errorAApprox = zeros(length(paramARR),length(kapARR));

kappa = 2;
sigma1 = 0.9;
sigma2 = 0.85;
mu = 0;
v0 = 6;
beta = 0.999;
syms omega D;
gamma = []; 
Omega = []; 
tic;
%% Find all solutions of polynomial, keep only unique and actual solutions

for j = 1:length(kapARR)
    kappa = kapARR(j);
    gamma = [];Omega = []; 
    for i = 1:length(paramARR)
        k = paramARR(i); %%% ---------- change ----------
    
        % Dielectric function with kappa=2 (using mu=0, need to add mu*k later)
        if kappa==2
        D = @(omega) 1+(-1).*k.^(-2).*(2.*beta.*k.^2.*(2.*omega.^2+k.^2.*sigma1.^2).^(-3).*( ...
              4.*omega.^4+12.*k.^2.*omega.^2.*sigma1.^2+(sqrt(-1)*(-8)).*2.^(1/2).* ...
              k.^3.*omega.*sigma1.^3+(-3).*k.^4.*sigma1.^4)+(1/2).*(2.*((-1)+beta).* ...
              k.^2.*(k.^2.*sigma2.^2+2.*(omega+(-1).*k.*v0).^2).^(-3).*(3.*k.^4.* ...
              sigma2.^4+(sqrt(-1)*8).*2.^(1/2).*k.^3.*sigma2.^3.*(omega+(-1).*k.*v0)+( ...
              -12).*k.^2.*sigma2.^2.*(omega+(-1).*k.*v0).^2+(-4).*(omega+(-1).*k.*v0) ...
              .^4)+2.*((-1)+beta).*k.^2.*(k.^2.*sigma2.^2+2.*(omega+k.*v0).^2).^(-3).* ...
              (3.*k.^4.*sigma2.^4+(sqrt(-1)*8).*2.^(1/2).*k.^3.*sigma2.^3.*(omega+k.* ...
              v0)+(-12).*k.^2.*sigma2.^2.*(omega+k.*v0).^2+(-4).*(omega+k.*v0).^4))); 
        elseif kappa==4
        D = @(omega) 1+(-1).*k.^(-2).*(2.*beta.*k.^2.*(2.*omega.^2+5.*k.^2.*sigma1.^2).^(-5) ...
              .*(16.*omega.^8+224.*k.^2.*omega.^6.*sigma1.^2+1400.*k.^4.*omega.^4.* ...
              sigma1.^4+7000.*k.^6.*omega.^2.*sigma1.^6+(sqrt(-1)*(-3200)).*10.^(1/2) ...
              .*k.^7.*omega.*sigma1.^7+(-4375).*k.^8.*sigma1.^8)+(1/2).*(2.*((-1)+ ...
              beta).*k.^2.*(5.*k.^2.*sigma2.^2+2.*(omega+(-1).*k.*v0).^2).^(-5).*( ...
              4375.*k.^8.*sigma2.^8+(-7000).*k.^6.*sigma2.^6.*(omega+(-1).*k.*v0).^2+( ...
              -1400).*k.^4.*sigma2.^4.*(omega+(-1).*k.*v0).^4+(-224).*k.^2.* ...
              sigma2.^2.*(omega+(-1).*k.*v0).^6+(-16).*(omega+(-1).*k.*v0).^8+(sqrt( ...
              -1)*(-3200)).*10.^(1/2).*k.^7.*sigma2.^7.*((-1).*omega+k.*v0))+2.*((-1)+ ...
              beta).*k.^2.*(5.*k.^2.*sigma2.^2+2.*(omega+k.*v0).^2).^(-5).*(4375.* ...
              k.^8.*sigma2.^8+(sqrt(-1)*3200).*10.^(1/2).*k.^7.*sigma2.^7.*(omega+k.* ...
              v0)+(-7000).*k.^6.*sigma2.^6.*(omega+k.*v0).^2+(-1400).*k.^4.* ...
              sigma2.^4.*(omega+k.*v0).^4+(-224).*k.^2.*sigma2.^2.*(omega+k.*v0).^6+( ...
              -16).*(omega+k.*v0).^8)));
        end
    
        Y = double(vpasolve(D(omega)==0,omega)) + mu*k;
    
        gamma = [gamma, imag(Y)];
        Omega = [Omega, real(Y)]; 
    end
    numSol = length(gamma(:,1));
    
    % Run to get only unique solutions
    [gamma, Omega, numSol] = getUnique(gamma, Omega);
    
    % Run to get only actual solutions
    % [gamma, Omega, numSol] = getActualSol(gamma, Omega, D, paramARR);
    
    %% Compute the analytic approximation and Xie/Weideman root, compute 
    % relative error between them and root from above with largest imaginary
    % part
    for i = 1:length(paramARR)
        k = paramARR(i);
    
        Y = Omega(:,i)+1i*gamma(:,i);
        if kappa ==2; temp1 = Y(imag(Y)<0.601044590483068*0.95*k-3.44069961804043e-07);
        elseif kappa ==4; temp1 = Y(imag(Y)<1.34278881453794*0.95*k+0.000160077999823189); end

        temp2 = sort(1i*temp1,'ComparisonMethod','real'); 
        exact = -1i*temp2(1); % choose root with largest imaginary part
        omegaExact(i,j) = abs(real(exact))+1i*imag(exact);

        init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k,sigma1,sigma2,0,v0,beta,kappa) + mu*k;

        % init_guess = omegaExact(i,j);
        % init_guess = BohmGross_KapBoT(k,sigma1,sigma2,v0,beta,kappa) + mu*k;
        xi_guess_shiftscaled = (init_guess/k-mu)/sigma1;
        omegaXieShiftScaled(i,j) = BiKappa_Disp_Using_Xie(k*sigma1,1,sigma2/sigma1,0,v0/sigma1,beta,kappa,xi_guess_shiftscaled)*sigma1*k + mu*k;
        omegaAApprox(i,j) = AnalyticApproxKappa(k,sigma1,sigma2,v0,beta,kappa) + mu*k;

        errorXie(i,j) = abs(imag(omegaExact(i,j))-imag(omegaXieShiftScaled(i,j)))./abs(imag(omegaExact(i,j)));
        errorAApprox(i,j) = abs(imag(omegaExact(i,j))-imag(omegaAApprox(i,j)))./abs(imag(omegaExact(i,j)));
    end
end
toc;

%% Run Figures
close all;
colorSet = {'#dc143c','#cb2888','#9932cc','#6c4ce2','#2f5ada'};%red->blue
colors = {'#dc143c','#dc143c','#2f5ada','#2f5ada'};

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% % imaginary part of roots
% figure
% plot(paramARR,gamma,'.-','linewidth',1); hold on
% plot(paramARR,imag(omegaAApprox(:,1)),'.k','linewidth',1,'MarkerSize',6);
% title("Kappa Bump-on-Tail: All Roots $\gamma(k)$",'FontSize',14)
% lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
% title(lgd,'Solution #'); 
% xlabel("$k$"); ylabel("$\gamma(k)$"); grid on;
% xlim([0.001,0.275]);  ylim([-0.3,0.2])
% set(gcf, 'PaperPosition', [0 0 6 5]); %Position the plot further to the upper-left conder
% set(gcf, 'PaperSize', [6 5]); % Extends the plot to fill the entire paper
% saveas(gcf, 'Figs/LevelCurves_Omegak_ALL.pdf')
% 
% figure
% plot(paramARR,Omega,'.-','linewidth',1); hold on
% plot(paramARR,real(omegaAApprox(:,1)),'.k','linewidth',1,'MarkerSize',6);
% title("Kappa Bump-on-Tail: All Roots $\Omega(k)$",'FontSize',14)
% lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
% title(lgd,'Solution #'); 
% xlabel("$k$"); ylabel("$\Omega(k)$"); grid on;
% xlim([0.001,0.275]); ylim([0,1.8])
% set(gcf, 'PaperPosition', [0 0 6 5]); %Position the plot further to the upper-left conder
% set(gcf, 'PaperSize', [6 5]); % Extends the plot to fill the entire paper
% saveas(gcf, 'Figs/LevelCurves_gammak_ALL.pdf')

figure
% imaginary part of roots
plot(paramARR,imag(omegaExact(:,1)),'-k','linewidth',1); hold on
plot(paramARR,imag(omegaXieShiftScaled(:,1)),'s','linewidth',1,'MarkerSize',5);
plot(paramARR,imag(omegaAApprox(:,1)),'.','linewidth',1,'MarkerSize',10);
plot(paramARR,imag(omegaExact(:,2)),'--k','linewidth',1);
plot(paramARR,imag(omegaXieShiftScaled(:,2)),'s','linewidth',1,'MarkerSize',5);
plot(paramARR,imag(omegaAApprox(:,2)),'.','linewidth',1,'MarkerSize',10);
title("{Kappa Bump-on-Tail Level Curves}",'FontSize',16)
    xlabel("$k$",'FontSize',14); ylabel("$\gamma(k)$",'FontSize',14);
    set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
    set(gca,'Yscale','linear','XScale','log')
    colororder(colorSet([1,1,1,1,5,5]))
    set(gcf, 'PaperPosition', [0 0 6 4.5]); %Position the plot further to the upper-left conder
    set(gcf, 'PaperSize', [6 4.5]); % Extends the plot to fill the entire paper
    % plot(0.9976/v0*[1 1],get(gca,'ylim'),'k--')
    % plot(1.0132/v0*[1 1],get(gca,'ylim'),'k--')
    % plot(paramARR,real(omegaExact),'k')
    xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.8,1])
    legend('Exact $(\kappa\!=\!2)$','Weideman $(\kappa\!=\!2)$','Sarkar $(\kappa\!=\!2)$',...
        'Exact $(\kappa\!=\!4)$','Weideman $(\kappa\!=\!4)$','Sarkar $(\kappa\!=\!4)$',...
        '$k\approx \Omega/v_0$','location','SouthWest','FontSize',12)
    saveas(gcf, 'Figs/gammak_Sarkar.pdf')
    saveas(gcf, 'Figs/gammak_Sarkar.fig')

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

% % relative error between largest root and Weideman/Xie and Analytic Approx
% figure
% loglog(paramARR,errorXie(:,1),'.','linewidth',1,'MarkerSize',10); hold on
% loglog(paramARR,errorAApprox(:,1),'--','linewidth',1); 
% loglog(paramARR,errorXie(:,2),'.','linewidth',1,'MarkerSize',10); 
% loglog(paramARR,errorAApprox(:,2),'--','linewidth',1); 
% title('Relative Error in Growth Rate','FontSize',16);
%     legend('Weideman $(\kappa\!=\!2)$','Analytic Approx. $(\kappa\!=\!2)$',...
%         'Weideman $(\kappa\!=\!4)$','Analytic Approx. $(\kappa\!=\!4)$','Location','Best','FontSize',12);
%     xlabel('$k$','FontSize',14);
%     ylabel('$\big|\frac{\gamma-\gamma_{i}}{\gamma}\big|$ (log scale)','FontSize',14);
%     set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
%     set(gca,'Yscale','log','XScale','log')
%     colororder(colors);
%     set(gcf, 'PaperPosition', [0 0 6 4.5]); %Position the plot further to the upper-left conder
%     set(gcf, 'PaperSize', [6 4.5]); % Extends the plot to fill the entire paper
%     saveas(gcf, 'Figs/gammaRelErrork_Sarkar.pdf')
%     saveas(gcf, 'Figs/gammaRelErrork_Sarkar.fig')

    save('Comparison_SarkarAnalytic_DBoTkappa_data.mat')


%% Functions for eliminating non-physical and duplicate roots
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
    % z44 = floor(colMax*0.44); % index at 44% of length_k
    % z42 = floor(colMax*0.42); % index at 42% of length_k
    for j=1:numSol
        Fzmax = max( double(abs(D(gamma(j,:)+1i*Omega(j,:)))) );
        dgammaEnd = (gamma(j,z90)-gamma(j,z88))/(paramARR(z90)-paramARR(z88)); %slope of gamma
        % dgammaMid = (gamma(j,z44)-gamma(j,z42))/(paramARR(z44)-paramARR(z42)); %slope of gamma
        if Fzmax < tol && dgammaEnd < dtol && dgammaEnd ~= 0.0 %&& dgammaMid < dtol*2
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