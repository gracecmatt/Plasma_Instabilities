%% Incomplete Bi-Maxwellian Test Script for Accuracy of Methods
% Runs level curves in each variable and generates N^2 random samples
% inside variation limits. Plots level curves and relative error in omega.
% This script uses the dispersion function of the Bi-Maxwellian as the "true" solution, the
% spectral code for the initial guess, and the Xie code shifted, not scaled.
clear; clc;

N = 100; % number of samples in level curves
Nparams = 7;
M = 8; % slope of erf approximation to step function, for spectral code

% parameters = [k, sigma1, sigma2, mu1, mu2, beta, nu]
baseparams = [0.5; 1; 0.6; 1; 4; 0.9; -20];
txtparams = ["k","\sigma_1","\sigma_2","\mu_1","\mu_2","\beta","\nu"];

% set variations
var = 0.5; % x 100% variation considered
xl = (1-var)*baseparams;
xu = (1+var)*baseparams;
xu(6) = min(xu(6),0.97); % keep beta < 1
xl(6) = max(xl(6),0.58); % keep beta > 0.5

% Pre-allocate memory
paramsarr = zeros(Nparams,N);
omega.spectral = zeros(Nparams,N);
omega.xie = zeros(Nparams,N);
omega.xie_shift = zeros(Nparams,N);
omega.xie_shiftscaled = zeros(Nparams,N);
omega.dielectric = zeros(Nparams,N);
errorRand.omega = zeros(1,N^2);
errorRand.spectral = zeros(1,N^2);

%% Run level curves in each parameter
% initialize parameters to base 
params = baseparams;

% iterate through each parameter, set its array, and compute level curves
for i = 1:Nparams
    % fix upper/lower bound if base parameter is zero
    if baseparams(i)==0
        xl(i) = -var;
        xu(i) = var;
    end
    % define parameter array around base using variation
    paramsarr(i,:) = linspace(xl(i),xu(i),N);

    % iterate ith parameter through its parameter array
    for j=1:N
        params(i) = paramsarr(i,j);
            k=params(1);
            sigma1=params(2);
            sigma2=params(3);
            mu1=params(4);
            mu2=params(5);
            beta=params(6);
            nu=params(7);

        init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma1,sigma2,0,mu2-mu1,beta,nu,M) + mu1*k;
        % init_guess = BohmGross_BiMax(k,sigma1,sigma2,0,mu2-mu1,beta,nu) + mu1*k;

        xi_guess = init_guess/k;
        xi_guess_shift = init_guess/k-mu1;
        xi_guess_shiftscaled = (init_guess/k-mu1)/sigma1;

        omega.spectral(i,j) = init_guess; 
        omega.xie(i,j) = IncompleteBiMax_Disp_Using_Xie(k,sigma1,sigma2,mu1,mu2,beta,nu,xi_guess)*k;
        omega.xie_shift(i,j) = IncompleteBiMax_Disp_Using_Xie(k,sigma1,sigma2,0,mu2-mu1,beta,nu-mu1,xi_guess_shift)*k + mu1*k;
        omega.xie_shiftscaled(i,j) = IncompleteBiMax_Disp_Using_Xie(k*sigma1,1,sigma2/sigma1,0,(mu2-mu1)/sigma1,beta,(nu-mu1)/sigma1,xi_guess_shiftscaled)*sigma1*k + mu1*k;
        omega.dielectric(i,j) = BiMax_dielectric([k,sigma1,sigma2,0,mu2-mu1,beta],xi_guess-mu1)*k + mu1*k;
    end

    % re-save params as base parameters for next level curve
    params = baseparams;
end

% calculate relative error
errorRe.spectral = abs(real(omega.dielectric)-real(omega.spectral))./abs(real(omega.dielectric));
errorRe.xie = abs(real(omega.dielectric)-real(omega.xie))./abs(real(omega.dielectric));
errorRe.xie_shift = abs(real(omega.dielectric)-real(omega.xie_shift))./abs(real(omega.dielectric));
errorRe.xie_shiftscaled = abs(real(omega.dielectric)-real(omega.xie_shiftscaled))./abs(real(omega.dielectric));

errorIm.spectral = abs(imag(omega.dielectric)-imag(omega.spectral))./abs(imag(omega.dielectric));
errorIm.xie = abs(imag(omega.dielectric)-imag(omega.xie))./abs(imag(omega.dielectric));
errorIm.xie_shift = abs(imag(omega.dielectric)-imag(omega.xie_shift))./abs(imag(omega.dielectric));
errorIm.xie_shiftscaled = abs(imag(omega.dielectric)-imag(omega.xie_shiftscaled))./abs(imag(omega.dielectric));

%% Plot Level Curves and Percent Errors
% close all
% txtbase = string(['$(k,\sigma_1,\sigma_2,\mu_1,\mu_2,\beta)$=(',num2str(baseparams(1)),...
%     ',',num2str(baseparams(2)),',',num2str(baseparams(3)),',',num2str(baseparams(4)),...
%     ',',num2str(baseparams(5)),',',num2str(baseparams(6)),',',num2str(baseparams(7)),')']);
% txtleg1 = {'Spectral','Xie','Shifted Xie','Shifted & Scaled Xie','Dielectric Func.'};
% txtleg2 = {'Spectral','Xie','Shifted Xie','Shifted & Scaled Xie'};
% 
% newcolors = {'#4363d8','#e6194B','#3cb44b','#7E2F8E','#42d4f4'};
% 
% % plot the level curves
% fig1 = figure; tiledlayout(2,4);
% for i=1:Nparams
%     nexttile
%     plot(paramsarr(i,:), imag(omega.spectral(i,:)),'.-'); hold on
%     plot(paramsarr(i,:), imag(omega.xie(i,:)),'.-');
%     plot(paramsarr(i,:), imag(omega.xie_shift(i,:)),'.-');
%     plot(paramsarr(i,:), imag(omega.xie_shiftscaled(i,:)),'.-');
%     plot(paramsarr(i,:), imag(omega.dielectric(i,:)),'k');
%         title('$\gamma('+txtparams(i)+')$ Level Curves','Interpreter','latex','FontSize',14);
%         xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
%         ylabel('$\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',12);
%         colororder(newcolors);
% end
%     leg = legend(txtleg1,'Orientation','Horizontal','FontSize',8);
%     leg.Layout.Tile = 'north';
%     txt = 'Incomplete BiMaxwellian $\gamma(p)$ around';
%     sgtitle({txt,txtbase},'Interpreter','latex','FontSize',16);
% 
% fig2 = figure; tiledlayout(2,4);
% for i=1:Nparams
%     nexttile
%     plot(paramsarr(i,:), real(omega.spectral(i,:)),'.-'); hold on
%     plot(paramsarr(i,:), real(omega.xie(i,:)),'.-');
%     plot(paramsarr(i,:), real(omega.xie_shift(i,:)),'.-');
%     plot(paramsarr(i,:), real(omega.xie_shiftscaled(i,:)),'.-');
%     plot(paramsarr(i,:), real(omega.dielectric(i,:)),'k');
%         title('$\Omega('+txtparams(i)+')$ Level Curves','Interpreter','latex','FontSize',14);
%         xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
%         ylabel('$\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',12);
%         colororder(newcolors);
% end
%     leg = legend(txtleg1,'Orientation','Horizontal','FontSize',8);
%     leg.Layout.Tile = 'north';
%     txt = 'Incomplete BiMaxwellian $\Omega(p)$ around';
%     sgtitle({txt,txtbase},'Interpreter','latex','FontSize',16);
% 
% % plot the percent error
% fig3 = figure; tiledlayout(2,4);
% for i=1:Nparams
%     nexttile
%     semilogy(paramsarr(i,:), errorIm.spectral(i,:),'.-'); hold on
%     semilogy(paramsarr(i,:), errorIm.xie(i,:),'.-');
%     semilogy(paramsarr(i,:), errorIm.xie_shift(i,:),'.-');
%     semilogy(paramsarr(i,:), errorIm.xie_shiftscaled(i,:),'.-');
%         title('Relative Error in $\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',14);
%         xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
%         ylabel('Relative Error','Interpreter','latex','FontSize',12);
%         ylim([0,1]); grid on; colororder(newcolors);
% end
%     leg = legend(txtleg2,'Orientation','Horizontal','FontSize',9);
%     leg.Layout.Tile = 'north';
%     txt = 'Incomplete BiMaxwellian Error in $\gamma(p)$ around';
%     sgtitle({txt,txtbase},'Interpreter','latex','FontSize',16);
% 
% fig4 = figure; tiledlayout(2,4);
% for i=1:Nparams
%     nexttile
%     semilogy(paramsarr(i,:), errorRe.spectral(i,:),'.-'); hold on
%     semilogy(paramsarr(i,:), errorRe.xie(i,:),'.-');
%     semilogy(paramsarr(i,:), errorRe.xie_shift(i,:),'.-');
%     semilogy(paramsarr(i,:), errorRe.xie_shiftscaled(i,:),'.-');
%         title('Relative Error in $\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',14);
%         xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
%         ylabel('Relative Error','Interpreter','latex','FontSize',12);
%         ylim([0,1]); grid on; colororder(newcolors);
% end
%     leg = legend(txtleg2,'Orientation','Horizontal','FontSize',9);
%     leg.Layout.Tile = 'north';
%     txt = 'Incomplete BiMaxwellian Error in $\Omega(p)$ around';
%     sgtitle({txt,txtbase},'Interpreter','latex','FontSize',16);
% 
% pause(0.5) % pause to let figures show up
%% Run simulations N^2 times with randomly selected parameters
rng('shuffle');
rng(sum(100*clock));
Xs = 2*rand(N^2,Nparams) - 1; % do sampling in serial

% tic;
% if N==50; parpool(5); elseif N==100; parpool(10); else; parpool(12); end
parpool(10);
parfor j = 1:N^2
    % Randomly sample parameters within acceptable ranges
    randparams = 1/2*(diag(xu - xl)*Xs(j,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with randomly drawn parameters
    % init_guess = BohmGross_BiMax(randparams(1),randparams(2),randparams(3),0,randparams(5)-randparams(4),randparams(6));
    init_guess = Vlasov_1D_linearized_Steve_v4(randparams(1),randparams(2),randparams(3),0,randparams(5)-randparams(4),randparams(6),randparams(7),M);
    xi_guess = init_guess/randparams(1); % shifted (not scaled)
 
    spectral = init_guess + randparams(4)*randparams(1);
    xie = IncompleteBiMax_Disp_Using_Xie(randparams(1),randparams(2),randparams(3),0,randparams(5)-randparams(4),randparams(6),randparams(7)-randparams(4),xi_guess)*randparams(1) + randparams(4)*randparams(1);
    dielectric = BiMax_dielectric([randparams(1),randparams(2),randparams(3),0,randparams(5)-randparams(4),randparams(6)],xi_guess)*randparams(1) + randparams(4)*randparams(1);
        
    errorRand(j).omega = abs(real(dielectric)-real(xie))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(xie))/abs(imag(dielectric));
    errorRand(j).spectral = abs(real(dielectric)-real(spectral))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(spectral))/abs(imag(dielectric));
end
delete(gcp('nocreate'));
errorMax_spectral = max(abs([errorRand.spectral]));
errorMax = max(abs([errorRand.omega]));
errorMean = mean(abs([errorRand.omega]));
% toc

%% Plot the percent error in the N^2 randomly drawn samples
% fig5 = figure; tiledlayout(1,2);
% nexttile
% semilogy(real([errorRand.omega]),'*');
%     title('$\Omega$ Relative Error','Interpreter','latex','FontSize',14);
%     xlabel('Sample \#','Interpreter','latex','FontSize',14);
%     ylabel('$\big|\frac{\Omega_{dielctric}-\Omega_{Xie}}{\Omega_{dielectric}}\big|$','Interpreter','latex','FontSize',14);
%     xlim([0,N^2]); grid on
% nexttile
% semilogy(imag([errorRand.omega]),'*');
%     title('$\gamma$ Relative Error','Interpreter','latex','FontSize',14);
%     xlabel('Sample \#','Interpreter','latex','FontSize',14);
%     ylabel('$\big|\frac{\gamma_{dielectric}-\gamma_{Xie}}{\gamma_{dielectric}}\big|$','Interpreter','latex','FontSize',14);
%     xlim([0,N^2]); grid on
% 
% txt1 = ['Incomplete BiMaxwellian, ',num2str(var*100),'\% variation on'];
% txt2 = ['$(k,\sigma_1,\sigma_2,\mu_1,\mu_2,\beta,\nu)$=(',...
%     num2str(baseparams(1)),',',num2str(baseparams(2)),',',num2str(baseparams(3)),',',num2str(baseparams(4)),',',num2str(baseparams(5)),',',num2str(baseparams(6)),',',num2str(baseparams(7)),')'];
% sgtitle({txt1;txt2},'Interpreter','latex','FontSize',16)
%     txt = ['Maximum Relative Error: ',num2str(errorMax,'%0.4e')];
%     leg = legend('','Orientation','Horizontal','color','none'); legend boxoff;
%     leg.Layout.Tile = 'north';
%     annotation('textbox',[.25 .65 .6 .2], ...
%         'String',txt,'EdgeColor','none','Interpreter','latex','FontSize',12);

%% Save the testing data and figures
% eps is a vector format file type that works well in LaTeX
% pdf is a vector format file type that is easy to share and view

save(['Data\TestingDataIBiMax_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '_data.mat'])

% savefig(fig1,  ['Figs\LevelCurvesBiMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig1, ['Figs\LevelCurvesBiMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig1, ['Figs\LevelCurvesBiMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% savefig(fig2,  ['Figs\LevelCurvesBiMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig2, ['Figs\LevelCurvesBiMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig2, ['Figs\LevelCurvesBiMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% savefig(fig3,  ['Figs\RelErrorBiMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig3, ['Figs\RelErrorBiMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig3, ['Figs\RelErrorBiMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% savefig(fig4,  ['Figs\RelErrorBiMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig4, ['Figs\RelErrorBiMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig4, ['Figs\RelErrorBiMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% savefig(fig5,  ['Figs\SampRandBiMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.fig'])
% hgexport(fig5, ['Figs\SampRandBiMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig5, ['Figs\SampRandBiMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% clear fig1 fig2 fig3 fig4 fig5
% save(['Data\TestingDataBiMax_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '_data.mat'])
