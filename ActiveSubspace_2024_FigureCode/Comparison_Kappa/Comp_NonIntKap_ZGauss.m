%% Non-Integer Kappa Comparison - Mace & Hellburg (1995) 
clear; clc;

% initialize kappas and other parameters/vectors
% assuming: mu = 0; sigma = 1; k = 0.5;
kapARR = round(logspace(0.18,0.7,20),2);
N = 1000; % number of samples
Zp = zeros(N,length(kapARR)); % results from Weideman algorithm
Zpgauss = zeros(N,length(kapARR)); % results from Gauss HyperGeometric function
error = zeros(N,length(kapARR));

% set z values to evaluate at from measurements of gamma from other simulations
zRe = linspace(-0.1,2,N)';
zIm = -1i*logspace(-8,-1,N)';
z = zRe(randsample(N,N),:) + zIm(randsample(N,N),:); % randomly combine real/imaginary parts

% iterate through all kappa values
for j = 1:length(kapARR)
    kappa = kapARR(j);
    
    % ========= compute Z, Zp using Fourier series approximation ==========
    F = ['(gamma(',num2str(kappa,16),')/gamma(',num2str(kappa,16),'-0.5)/sqrt(pi*',...
        num2str(kappa,16),'))*(1+v.^2/',num2str(kappa,16),').^(-',num2str(kappa,16),'-1)'];
    Fn = 0; % 0 = option to define your own F
    Nfourier = 1000; % number of Fourier coefficients to take
    [Zp(:,j),~] = zetaph(z, Fn, F, Nfourier); % dielectric function
    
    % ========= compute Z, Zp using Gauss hypergeometric function ==========
    x = 1/2*(1-z/(1i*sqrt(kappa)));
    % Zgauss = 1i*(kappa+1/2)*(kappa-1/2)/(kappa^(3/2)*(kappa+1))*hypergeom( [1, 2*kappa+2], kappa+2, x );
    Zpgauss(:,j) = -(kappa+1/2)*(kappa-1/2)/(kappa^2*(kappa+2))*hypergeom( [2, 2*kappa+3], kappa+3, x );

    error(:,j) = abs(Zp(:,j) - Zpgauss(:,j))./abs(Zpgauss(:,j));
end
%% 

Zmax = max(error,[],'all');
Zmean = mean(error,'all');
Zstd = std(error,0,'all');

%% ========= plot =========
close all;
color = {'#dc143c','#dc3798','#a700f7','#6344fc','#5684ff','#2f5ada'};%red->blue
color2 = {'#f21642','#e71a4b','#dc1d54','#d2205d','#c82465','#bd276e','#b32a76','#a92d7f','#9e3187','#943490','#7e3ca3','#733fac','#6943b5','#5f46bd','#5549c6','#4a4cce','#404fd7','#3653df','#2247df','#1548e8','#084af0'};
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

semilogy(kapARR,mean(error),'k-','linewidth',1); hold on
for j = 1:length(kapARR)
    semilogy(kapARR(j),mean(error(:,j)),'*','linewidth',2); 
    grid on; ylim([4*10^(-16),10^(-14)]); xlim([1.45,5.15])
end
    colororder(color2);
    xlabel("Kappa Value $(\kappa)$",'FontSize',14);
    ylabel("Mean Relative Error $\left( \frac{|Z_p-Z_{Gauss}|}{|Z_{Gauss}|} \right)$",'FontSize',14)
    title("Mean Relative Error in Weideman Algorithm vs. $\kappa$",'FontSize',16)
    set(gcf, 'PaperPosition', [0 0 6 4.5]); %Position the plot further to the upper-left conder
    set(gcf, 'PaperSize', [6 4.5]); % Extends the plot to fill the entire paper
    saveas(gcf, 'Figs/ZpCompare_smallKapLimit.pdf')
    saveas(gcf, 'Figs/ZpCompare_smallKapLimit.fig')


%%
loglog(1i*zIm*0.5,abs(Zp-Zpgauss),'linewidth',1.5,'markersize',5);
title('Relative Error in Modified Plasma Z Function','FontSize',16)
    ylabel('$\frac{\left|Z_p-Z_{p,W} \right|}{|Z_p|}$','FontSize',14)
    xlabel('Sample Values of $-\gamma$','FontSize',14)
    legend("$\kappa\!=\!"+num2str(kapARR(1))+"$","$\kappa\!=\!"+num2str(kapARR(2))+"$",...
        "$\kappa\!=\!"+num2str(kapARR(3))+"$","$\kappa\!=\!"+num2str(kapARR(4))+"$",...
        "$\kappa\!=\!"+num2str(kapARR(5))+"$","$\kappa\!=\!"+num2str(kapARR(6))+"$",'location','southeast')
    set(gca, 'YGrid', 'on', 'XGrid', 'off', 'YMinorGrid', 'off');
    set(gcf, 'PaperPosition', [0 0 6 4.5]); %Position the plot further to the upper-left conder
    set(gcf, 'PaperSize', [6 4.5]); % Extends the plot to fill the entire paper
    xlim([min(1i*zIm*0.5),max(1i*zIm*0.5)]);
    ax = gca;
    colororder(color);
    ax.LineStyleOrder = ["o";"+";"*";"s";"d";"h";">";"v"];
    ax.LineStyleCyclingMethod = "withcolor";
    saveas(gcf, 'Figs/ZpCompare_nonIntKap.pdf')
    saveas(gcf, 'Figs/ZpCompare_nonIntKap.fig')