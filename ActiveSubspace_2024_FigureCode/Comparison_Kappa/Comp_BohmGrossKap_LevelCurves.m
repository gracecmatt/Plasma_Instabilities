%% Kappa Distribution Comparison Level Curves
% GM 1/11/24

clear; clc; 
%%% --------- change these values with parameter for level curve ---------
% paramARR = linspace(0.005,0.75,150);
paramARR = 0.5*logspace(-2.3,0,500);
% paramARR=paramARR(10:end);
xaxis = "k"; 
filenameEnding = "k.png"; 
%%% ---------- end change section ----------
kapARR = [2,4];
% initialize vectors
omegaExact = zeros(length(kapARR),length(paramARR));
omegaXieShiftScaled = zeros(length(kapARR),length(paramARR));
omegaBohmGross = zeros(length(kapARR),length(paramARR));
errorXie = zeros(length(kapARR),length(paramARR));
errorBG = zeros(length(kapARR),length(paramARR));

% choose parameters
k = 0.5;
sigma = 1.5;
mu = 1;

% initialize for symbolic root-finder vpasolve()
syms omega D;
gamma = []; 
Omega = []; 

%% find the "actual" solutions -> polynomial roots
for j = 1:length(kapARR)
    kappa = kapARR(j);
    gamma = []; Omega = [];
    for i = 1:length(paramARR)
        k=paramARR(i); %%% ---------- change level curve parameter ----------
    
        % dielectric functions exported from Wolfram Mathematica (w/ zero mean)
        if (kappa==2)
            D(omega) = 1+(-2).*k.^(-2).*(2.*omega.^2+k.^2.*sigma.^2).^(-3).*(4.*k.^2.*omega.^4+ ...
              12.*k.^4.*omega.^2.*sigma.^2+(sqrt(-1)*(-8)).*2.^(1/2).*k.^5.*omega.* ...
              sigma.^3+(-3).*k.^6.*sigma.^4);
        elseif (kappa==4)
            D(omega) = 1+(-2).*k.^(-2).*(2.*omega.^2+5.*k.^2.*sigma.^2).^(-5).*(16.*k.^2.* ...
              omega.^8+224.*k.^4.*omega.^6.*sigma.^2+1400.*k.^6.*omega.^4.*sigma.^4+ ...
              7000.*k.^8.*omega.^2.*sigma.^6+(sqrt(-1)*(-3200)).*10.^(1/2).*k.^9.* ...
              omega.*sigma.^7+(-4375).*k.^10.*sigma.^8);
        else
            % error statement here
        end
    
        % symbolic root finder finds all roots, Y is a complex vector
        % mu*k is added, applying invariance of dielectric function
        Y = double(vpasolve(D(omega)==0,omega)) + mu*k;
    
        gamma = [gamma, imag(Y)];
        Omega = [Omega, real(Y)];
    end
    numSol = length(gamma(:,1));
    
    % Run to get only unique solutions
    [gamma, Omega, numSol] = getUnique(gamma, Omega);
    
    % Run to get only actual solutions (real part is not zero or increasing for large k)
    [gamma, Omega, numSol] = getActualSol(gamma, Omega, D, paramARR);
    
    
    
    
    %% compute the solution using Xie/Weideman algorithm and Bohm-Gross approximation
    % compare to the polynomial root (above) with the largest imaginary part
    
    for i = 1:length(paramARR)
        k = paramARR(i);
    
        % find the root with the largest imaginary part
        Y = Omega(:,i)+1i*gamma(:,i);
        temp = sort(1i*Y,'ComparisonMethod','real');
        omegaExact(j,i) = -1i*temp(1);
    
        % Xie/Weideman algorithm
        init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k,sigma,0,kappa) + mu*k;
        xi_guess_shiftscaled = (init_guess/k-mu)/sigma; % shifted & scaled
        omegaXieShiftScaled(j,i) = Kappa_Disp_Using_Xie(k*sigma,1,0,kappa,xi_guess_shiftscaled)*sigma*k + mu*k;
    
        % Bohm-Gross approximation
        omegaBohmGross(j,i) = BohmGross_Kap(k,sigma,0,kappa) + mu*k;
        
        % compute the relative error in the imaginary parts
        errorXie(j,i) = abs(imag(omegaExact(j,i))-imag(omegaXieShiftScaled(j,i)));%./abs(imag(omegaExact(j,i)));
        errorBG(j,i) = abs(imag(omegaExact(j,i))-imag(omegaBohmGross(j,i)));%./abs(imag(omegaExact(j,i)));
    end

end
%%


Xiemax = max(errorXie,[],'all');
Xiemean = mean(errorXie,'all');
Xiestd = std(errorXie,0,'all');


BGmax = max(errorBG,[],'all');
BGmean = mean(errorBG,'all');
BGstd = std(errorBG,0,'all');

save('Data/Comparison_BohmGross_kappa_2_4_data.mat')


%% plot figures 
% with ideal formatting (as of 1/27/24)
close all;
% colorSet = {'#dc143c','#cb2888','#9932cc','#6c4ce2','#2f5ada'};%red->blue
colors = {'#000000','#dc143c','#dc143c','#000000','#2f5ada','#2f5ada'};
colors2 = {'#dc143c','#dc143c','#2f5ada','#2f5ada'};

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% plot the imaginary part of the solutions found using Xie algorithm and
% Bohm-Gross approximation
yExact = imag(omegaExact);
yXie = imag(omegaXieShiftScaled);
yBG = imag(omegaBohmGross);

% figure
% loglog(paramARR,yExact(1,:),'-','linewidth',1); hold on
% loglog(paramARR,yXie(1,:),'.','linewidth',1,'MarkerSize',8);
% loglog(paramARR,yBG(1,:),'s','linewidth',1,'MarkerSize',5)
% loglog(paramARR,yExact(2,:),'-','linewidth',1); 
% loglog(paramARR,yXie(2,:),'.','linewidth',1,'MarkerSize',8);
% loglog(paramARR,yBG(2,:),'s','linewidth',1,'MarkerSize',5);
% title('\textbf{Kappa $\gamma(k)$ Level Curves}','FontSize',16);
%     legend('Exact $(\kappa\!=\!2)$','Weideman $(\kappa\!=\!2)$',...
%         'Bohm-Gross $(\kappa\!=\!2)$','Exact $(\kappa\!=\!4)$',...
%         'Weideman $(\kappa\!=\!4)$','Bohm-Gross $(\kappa\!=\!4)$',...
%         'Location','NorthEast','FontSize',12);
%     xlabel('$k$','FontSize',14); xlim([paramARR(1),paramARR(end)]);
%     ylabel('$\gamma(k)$ (negative log)','FontSize',14);
%     set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
%     set(gca,'Yscale','log','XScale','Linear')
%     colororder(colors)
% 
%     set(gcf, 'PaperPosition', [0 0 6 4.5]); %Position the plot further to the upper-left conder
%     set(gcf, 'PaperSize', [6 4.5]); % Extends the plot to fill the entire paper
%     saveas(gcf, 'Figs/gammak_LevelCurves_BGWeideman.pdf')
%     saveas(gcf, 'Figs/gammak_LevelCurves_BGWeideman.fig')

% plot the relative error between Xie/Weideman solution, Bohm-Gross
% solution, and the polynomial root with the largest imaginary part
figure
plot(paramARR,log10(errorXie(1,:)),'s','linewidth',1,'MarkerSize',5); hold on
plot(paramARR,log10(errorBG(1,:)),'.','linewidth',1,'MarkerSize',10)
plot(paramARR,log10(errorXie(2,:)),'s','linewidth',1,'MarkerSize',5);
plot(paramARR,log10(errorBG(2,:)),'.','linewidth',1,'MarkerSize',10)
title('\textbf{Kappa $\gamma(k)$ Relative Error}','FontSize',16);
    legend('Weideman $(\kappa\!=\!2)$','Bohm-Gross $(\kappa\!=\!2)$',...
        'Weideman $(\kappa\!=\!4)$','Bohm-Gross $(\kappa\!=\!4)$',...
        'Location','Best','FontSize',12);
    xlabel('$k$','FontSize',14);
    ylabel('$\big|\frac{\gamma-\gamma_{i}}{\gamma}\big|$ (log scale)','FontSize',14);
    set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
    colororder(colors2)
    xlim([min(paramARR),max(paramARR)])
    set(gcf, 'PaperPosition', [0 0 6 4.5]); %Position the plot further to the upper-left conder
    set(gcf, 'PaperSize', [6 4.5]); % Extends the plot to fill the entire paper
    saveas(gcf, 'Figs/gammaRelErrork_Kap24_BGXie.pdf')
    saveas(gcf, 'Figs/gammaRelErrork_Kap24_BGXie.fig')

save('Data/Comparison_BohmGross_kappa_data.mat')



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