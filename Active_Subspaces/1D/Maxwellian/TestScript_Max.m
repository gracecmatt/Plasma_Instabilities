%% Maxwellian Test Script for Accuracy of Methods
% Runs level curves in each variable and generates N^2 random samples
% inside variation limits. Plots level curves and relative error in omega.
% This script uses the dispersion function as the "true" solution, the
% Bohm-Gross relation for the initial guess, and the Xie code shifted, not scaled.
clear; clc;

N = 100; % number of samples in level curves
Nparams = 3;

% parameters = [k, sigma, mu]
baseparams = [0.5; 2; 1];
txtparams = ["k","\sigma","\mu"];

% set variations
var = 0.50; % x 100% variation considered
xl = (1-var)*baseparams;
xu = (1+var)*baseparams;

% Pre-allocate memory
paramsarr = zeros(Nparams,N);
omega.spectral = zeros(Nparams,N);
omega.xie = zeros(Nparams,N);
omega.xie_shift = zeros(Nparams,N);
omega.xie_shiftscaled = zeros(Nparams,N);
omega.dielectric = zeros(Nparams,N);
Xs = zeros(N,Nparams);
errorRand.omega = zeros(1,N^2);
errorRand.spectral = zeros(1,N^2);

%% Run level curves in each parameter
% initialize parameters to base 
params = baseparams;

for i = 1:Nparams 
    % fix upper/lower bound if base parameter is zero
    if baseparams(i)==0
        xl(i) = -var;
        xu(i) = var;
    end
    % define parameter array around base using variation
    paramsarr(i,:) = linspace(xl(i),xu(i),N);

    for j=1:N 
        % iterate ith parameter through its parameter array
        params(i) = paramsarr(i,j);
            k=params(1);
            sigma=params(2);
            mu=params(3);
        
        init_guess = BohmGross_Max(k,sigma,0) + mu*k;
        % init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma,0) + mu*k;

        xi_guess = init_guess/k;
        xi_guess_shift = init_guess/k-mu;
        xi_guess_shiftscaled = (init_guess/k-mu)/sigma;

        omega.spectral(i,j) = init_guess;
        omega.xie(i,j) = Maxwellian_Disp_Using_Xie(k,sigma,mu,xi_guess)*k;
        omega.xie_shift(i,j) = Maxwellian_Disp_Using_Xie(k,sigma,0,xi_guess_shift)*k + mu*k;
        omega.xie_shiftscaled(i,j) = Maxwellian_Disp_Using_Xie(k*sigma,1,0,xi_guess_shiftscaled)*sigma*k + mu*k;
        omega.dielectric(i,j) = Max_dielectric([k,sigma,0],xi_guess-mu)*k + mu*k;
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
close all
txtbase = string(['$(k,\sigma,\mu)$=(',num2str(baseparams(1)),...
    ',',num2str(baseparams(2)),',',num2str(baseparams(3)),')']);
txtleg1 = {'Bohm-Gross','Xie','Shifted Xie','Shifted & Scaled Xie','Dielectric Func.'};
txtleg2 = {'Bohm-Gross','Xie','Shifted Xie','Shifted & Scaled Xie'};

newcolors = {'#4363d8','#e6194B','#3cb44b','#7E2F8E','#42d4f4'};

% plot the level curves
fig1 = figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    plot(paramsarr(i,:), imag(omega.spectral(i,:)),'.-'); hold on
    plot(paramsarr(i,:), imag(omega.xie(i,:)),'.-');
    plot(paramsarr(i,:), imag(omega.xie_shift(i,:)),'.-');
    plot(paramsarr(i,:), imag(omega.xie_shiftscaled(i,:)),'.-');
    plot(paramsarr(i,:), imag(omega.dielectric(i,:)),'k');
        title('$\gamma('+txtparams(i)+')$ Level Curves','Interpreter','latex','FontSize',14);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
        ylabel('$\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',12);
        colororder(newcolors);
end
    leg = legend(txtleg1,'Orientation','Horizontal','FontSize',8);
    leg.Layout.Tile = 'north';
    sgtitle('Maxwellian $\gamma(p)$ around '+txtbase,'Interpreter','latex','FontSize',16);

fig2 = figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    plot(paramsarr(i,:), real(omega.spectral(i,:)),'.-'); hold on
    plot(paramsarr(i,:), real(omega.xie(i,:)),'.-');
    plot(paramsarr(i,:), real(omega.xie_shift(i,:)),'.-');
    plot(paramsarr(i,:), real(omega.xie_shiftscaled(i,:)),'.-');
    plot(paramsarr(i,:), real(omega.dielectric(i,:)),'k');
        title('$\Omega('+txtparams(i)+')$ Level Curves','Interpreter','latex','FontSize',14);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
        ylabel('$\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',12);
        colororder(newcolors);
end
    leg = legend(txtleg1,'Orientation','Horizontal','FontSize',8);
    leg.Layout.Tile = 'north';
    sgtitle('Maxwellian $\Omega(p)$ around '+txtbase,'Interpreter','latex','FontSize',16);

% plot the relative error
fig3 = figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    semilogy(paramsarr(i,:), errorIm.spectral(i,:),'.-'); hold on
    semilogy(paramsarr(i,:), errorIm.xie(i,:),'.-');
    semilogy(paramsarr(i,:), errorIm.xie_shift(i,:),'.-');
    semilogy(paramsarr(i,:), errorIm.xie_shiftscaled(i,:),'.-');
        title('Relative Error in $\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',14);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
        ylabel('Relative Error','Interpreter','latex','FontSize',12);
        ylim([0,max(1,max(errorIm.spectral(i,:)))]); grid on; colororder(newcolors);
end
    leg = legend(txtleg2,'Orientation','Horizontal','FontSize',9);
    leg.Layout.Tile = 'north';
    sgtitle('Maxwellian Error in $\gamma(p)$ around '+txtbase,'Interpreter','latex','FontSize',16);

fig4 = figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    semilogy(paramsarr(i,:), errorRe.spectral(i,:),'.-'); hold on
    semilogy(paramsarr(i,:), errorRe.xie(i,:),'.-');
    semilogy(paramsarr(i,:), errorRe.xie_shift(i,:),'.-');
    semilogy(paramsarr(i,:), errorRe.xie_shiftscaled(i,:),'.-');
        title('Relative Error in $\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',14);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',12);
        ylabel('Relative Error','Interpreter','latex','FontSize',12);
        ylim([0,max(1,max(errorRe.spectral(i,:)))]); grid on; colororder(newcolors);
end
    leg = legend(txtleg2,'Orientation','Horizontal','FontSize',9);
    leg.Layout.Tile = 'north';
    sgtitle('Maxwellian Error in $\Omega(p)$ around '+txtbase,'Interpreter','latex','FontSize',16);

pause(0.5) % pause to let figures show up
%% Run simulations N^2 times with randomly selected parameters
rng('shuffle');
rng(sum(100*clock));
Xs = 2*rand(N^2,Nparams) - 1; % do sampling in serial

% tic;
% if N==50; parpool(5); elseif N==100; parpool(10); else; parpool(12); end % change this line
parpool(10);
parfor j = 1:N^2
    % Randomly sample parameters within acceptable ranges
    randparams = 1/2*(diag(xu - xl)*Xs(j,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with randomly drawn parameters
    init_guess = BohmGross_Max(randparams(1),randparams(2),0);
    xi_guess = init_guess/randparams(1); % shifted (not scaled)

    spectral = init_guess + randparams(3)*randparams(1);
    xie = Maxwellian_Disp_Using_Xie(randparams(1),randparams(2),0,xi_guess)*randparams(1) + randparams(3)*randparams(1);
    dielectric = Max_dielectric([randparams(1),randparams(2),0],xi_guess)*randparams(1) + randparams(3)*randparams(1);
    
    errorRand(j).omega = abs(real(dielectric)-real(xie))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(xie))/abs(imag(dielectric));
    errorRand(j).spectral = abs(real(dielectric)-real(spectral))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(spectral))/abs(imag(dielectric));
end
delete(gcp('nocreate'));
errorMax_spectral = max(abs([errorRand.spectral]));
errorMax = max(abs([errorRand.omega]));
errorMean = meanx(abs([errorRand.omega]));
% toc

%% Plot the percent error in the N^2 randomly drawn samples
fig5 = figure; tiledlayout(1,2);
nexttile
semilogy(real([errorRand.omega]),'*');
    title('$\Omega$ Relative Error','Interpreter','latex','FontSize',14);
    xlabel('Sample \#','Interpreter','latex','FontSize',14);
    ylabel('$\big|\frac{\Omega_{dielctric}-\Omega_{Xie}}{\Omega_{dielectric}}\big|$',...
        'Interpreter','latex','FontSize',14);
    xlim([0,N^2]); grid on
nexttile
semilogy(imag([errorRand.omega]),'*');
    title('$\gamma$ Relative Error','Interpreter','latex','FontSize',14);
    xlabel('Sample \#','Interpreter','latex','FontSize',14);
    ylabel('$\big|\frac{\gamma_{dielectric}-\gamma_{Xie}}{\gamma_{dielectric}}\big|$'...
        ,'Interpreter','latex','FontSize',14);
    xlim([0,N^2]); grid on

txt = ['Maxwellian, ',num2str(var*100),'\% variation on $(k,\sigma,\mu)$=(',...
    num2str(baseparams(1)),',',num2str(baseparams(2)),',',num2str(baseparams(3)),')'];
sgtitle(txt,'Interpreter','latex','FontSize',16)
    txt = ['Maximum Relative Error: ',num2str(errorMax,'%0.4e')];
    leg = legend('','Orientation','Horizontal','color','none'); legend boxoff;
    leg.Layout.Tile = 'north';
    annotation('textbox',[.25 .72 .6 .2], ...
        'String',txt,'EdgeColor','none','Interpreter','latex','FontSize',12)

%% Save the testing data and figures
% eps is a vector format file type that works well in LaTeX
% pdf is a vector format file type that is easy to share and view

save(['TestingDataMax_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '_data.mat'])

% savefig(fig1,  ['Figs\LevelCurvesMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig1, ['Figs\LevelCurvesMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig1, ['Figs\LevelCurvesMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% savefig(fig2,  ['Figs\LevelCurvesMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig2, ['Figs\LevelCurvesMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig2, ['Figs\LevelCurvesMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% savefig(fig3,  ['Figs\RelErrorMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig3, ['Figs\RelErrorMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig3, ['Figs\RelErrorMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% savefig(fig4,  ['Figs\RelErrorMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig4, ['Figs\RelErrorMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig4, ['Figs\RelErrorMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% savefig(fig5,  ['Figs\SampRandMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.fig'])
% hgexport(fig5, ['Figs\SampRandMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig5, ['Figs\SampRandMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% clear fig1 fig2 fig3 fig4 fig5
% save(['Data\TestingDataMax_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '_data.mat'])
