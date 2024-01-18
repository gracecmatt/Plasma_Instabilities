%% Kappa Distribution Comparison Level Curves
% GM 1/11/24

clear; clc; 
%%% --------- change these values with parameter for level curve ---------
paramARR = linspace(0.01,1,80);
xaxis = "k"; 
filenameEnding = "k.png"; 
%%% ---------- end change section ----------
% initialize vectors
omegaXieShiftScaled = zeros(1,length(paramARR));
omegaBohmGross = zeros(1,length(paramARR));
errorXie = zeros(1,length(paramARR));
errorBG = zeros(1,length(paramARR));

% choose parameters
k = 0.5;
sigma = 1;
mu = 1;
kappa = 2;

% initialize for symbolic root-finder vpasolve()
syms omega D;
gamma = []; 
Omega = []; 

%% find the "actual" solutions -> polynomial roots
for i = 1:length(paramARR)
    k=paramARR(i); %%% ---------- change level curve parameter ----------

    % dielectric functions exported from Wolfram Mathematica (w/ zero mean)
    if (kappa ==2)
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
    exact = -1i*temp(1);

    % Xie/Weideman algorithm
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k,sigma,0,kappa) + mu*k;
    xi_guess_shiftscaled = (init_guess/k-mu)/sigma; % shifted & scaled
    omegaXieShiftScaled(i) = Kappa_Disp_Using_Xie(k*sigma,1,0,kappa,xi_guess_shiftscaled)*sigma*k + mu*k;

    % Bohm-Gross approximation
    omegaBohmGross(i) = BohmGross_Kap(k,sigma,0,kappa) + mu*k;
    
    % compute the relative error in the imaginary parts
    errorXie(i) = abs(imag(exact)-imag(omegaXieShiftScaled(i)))./abs(imag(exact));
    errorBG(i) = abs(imag(exact)-imag(omegaBohmGross(i)))./abs(imag(exact));
end




%% plot figures 
% with ideal formatting (as of 1/11/24)
% close all;
newcolors = {'#4363d8','#e6194B','#3cb44b','#42d4f4','#f032e6','#000075',...
    '#fabed4','#f58231','#469990','#9A6324','#800000','#aaffc3','#a9a9a9'};

% plot the imaginary part of the frequency, Im(omega)=gamma
figure
plot(paramARR,gamma,'.-','linewidth',1)
title("Kappa: All Roots \gamma("+xaxis+")",'FontSize',14)
lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
title(lgd,'Solution #'); 
xlabel(xaxis); ylabel("\gamma("+xaxis+")");
grid on; colororder(newcolors)
% saveas(gcf,"Matlab/AllSolImages/gammaVs"+filenameEnding)
    
    % figure
    % plot(paramARR,Omega,'.-','linewidth',1)
    % title("Kappa: All Roots \Omega("+xaxis+")",'FontSize',14)
    % lgd = legend(strcat(num2str((1:numSol)')),'Location','eastoutside');
    % title(lgd,'Solution #');
    % xlabel(xaxis); ylabel("\Omega("+xaxis+")");
    % grid on; colororder(newcolors)
    % % saveas(gcf,"Matlab/AllSolImages/OmegaVs"+filenameEnding)

% plot the imaginary part of the solutions found using Xie algorithm and
% Bohm-Gross approximation
figure
plot(paramARR,imag(omegaXieShiftScaled),'.-','linewidth',1); hold on
plot(paramARR,imag(omegaBohmGross),'.-','linewidth',1)
title('\textbf{Kappa Distribution $\mathbf{\gamma(k)}$}','Interpreter','latex','FontSize',16);
    legend('Weideman/Xie Algorithm','Bohm-Gross Approximation','Location','SouthWest','Interpreter','latex','FontSize',12);
    xlabel('$k$','Interpreter','latex','FontSize',14);
    ylabel('$\gamma(k)$ (log scale)','Interpreter','latex','FontSize',14);
    grid on; colororder(newcolors)
    xlim([0.01,1])
    set(gcf, 'PaperPosition', [0 0 6 5]); %Position the plot further to the upper-left conder
    set(gcf, 'PaperSize', [6 5]); % Extends the plot to fill the entire paper
    % saveas(gcf, 'Figs/gammak_Kap4_BGXie.pdf')
    % saveas(gcf, 'Figs/gammak_Kap4_BGXie.fig')

% plot the relative error between Xie/Weideman solution, Bohm-Gross
% solution, and the polynomial root with the largest imaginary part
figure
plot(paramARR,log10(errorXie),'.-','linewidth',1); hold on
plot(paramARR,log10(errorBG),'.-','linewidth',1)
title('\textbf{Kappa $\mathbf{\gamma(k)}$ Relative Error}','Interpreter','latex','FontSize',16);
    legend('Weideman/Xie Algorithm','Bohm-Gross Approximation','Location','SouthWest','Interpreter','latex','FontSize',12);
    xlabel('$k$','Interpreter','latex','FontSize',14);
    ylabel('$\big|\frac{\gamma-\gamma_{i}}{\gamma}\big|$ (log scale)','Interpreter','latex','FontSize',14);
    grid on; colororder(newcolors)
    xlim([0.01,1])
    set(gcf, 'PaperPosition', [0 0 6 5]); %Position the plot further to the upper-left conder
    set(gcf, 'PaperSize', [6 5]); % Extends the plot to fill the entire paper
    % saveas(gcf, 'Figs/gammaRelErrork_Kap4_BGXie.pdf')
    % saveas(gcf, 'Figs/gammaRelErrork_Kap4_BGXie.fig')




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