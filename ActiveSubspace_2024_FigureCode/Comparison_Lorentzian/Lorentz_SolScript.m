%% 1D Lorentzian Solution Script Computing Relative Error in gamma
clear; clc;
%%% --- parameters to change when tuning the numerical algorithm ---
% L (zetaph) = 1
% M (spectral) = 2^10
% Vmax (spectral) = 140

N = 600;
paramARR(1,:) = linspace(0.2,1,N); % k array
paramARR(2,:) = linspace(0.5,4,N); % sigma array

errorXieIM = zeros(2,N); % relative error of imaginary part of root

mu = 10; % mean velocity

for i = 1:N % iterate through k array, keeping sigma constant
    sigma = 2;
    k = paramARR(1,i);

    exact = mu*k+1 + 1i.*(-sigma*k);

    init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma,0) + mu*k;
    xi_guess_shiftscaled = (init_guess/k-mu)/sigma;
    omegaXieShiftScaled =  Lorentzian_Disp_Using_Xie(k*sigma,1,0,xi_guess_shiftscaled)*sigma*k + mu*k;

    errorXieIM(1,i) = abs(imag(exact)-imag(omegaXieShiftScaled))./abs(imag(exact));
    if errorXieIM(2,i) == 0; errorXieIM(2,i) = 10^(-16); end
end
for i = 1:N % iterate through sigma array, keeping k constant
    k = 0.5;
    sigma = paramARR(2,i);

    exact = mu*k+1 + 1i.*(-sigma*k);

    init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma,0) + mu*k;
    xi_guess_shiftscaled = (init_guess/k-mu)/sigma;
    omegaXieShiftScaled =  Lorentzian_Disp_Using_Xie(k*sigma,1,0,xi_guess_shiftscaled)*sigma*k + mu*k;

    errorXieIM(2,i) = abs(imag(exact)-imag(omegaXieShiftScaled))./abs(imag(exact));
    if errorXieIM(2,i) == 0; errorXieIM(2,i) = 10^(-16); end
end

%% Run Relative Error Figure
close all;
y1 = log10(errorXieIM(1,:)); % blue: '#2f5ada'
y2 = log10(errorXieIM(2,:)); % red: '#dc143c'
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure; tiledlayout(1,1);
ax1 = axes('Position',[0.11 0.11 0.8 0.7]);
pl.errork = plot(ax1,paramARR(1,:),y1,'*','Color','#2f5ada','linewidth',1); 
    ax1.XColor = '#2f5ada'; ax1.Box = 'off'; hold off; 
    xlabel(ax1,'$k$','Interpreter','latex','FontSize',14);
    ylabel('$\big|\frac{\gamma-\gamma_{i}}{\gamma}\big|$ (log scale)','FontSize',14);
ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','color','none'); hold(ax2,'on')
pl.errorsig = plot(ax2,paramARR(2,:),y2,'*','Color','#dc143c','linewidth',1); 
    ax2.XColor = '#dc143c'; ax2.YColor = 'none'; ax2.Box = 'off'; hold off;
    xlabel(ax2,'$\sigma$','FontSize',14);
linkprop([ax1, ax2], {'ylim', 'Position'}); % Link the y limits and position together
    ylim([-14.9,-4.5]); grid on;
title('\textbf{Lorentzian $\gamma(k,\sigma)$ Relative Error}','FontSize',16);
    line(paramARR(1,1),0,'Parent',ax2,'Color','k'); % use to keep title visible

legend([pl.errorsig,pl.errork],{'Error in $\gamma(\sigma)$','Error in $\gamma(k)$'},'Location','SouthWest','FontSize',12);
    set(gcf, 'PaperPosition', [0 0 6 5]); %Position the plot further to the upper-left conder
    set(gcf, 'PaperSize', [6 5]); % Extends the plot to fill the entire paper
    saveas(gcf, "Figs/gammaRelError_Lorentz.pdf")
    saveas(gcf, "Figs/gammaRelError_LorentzFIG.fig")