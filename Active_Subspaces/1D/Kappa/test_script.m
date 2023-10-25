%% Kappa Test Script for Accuracy of Methods
% Runs level curves in each variable and generates N^2 random samples
% inside variation limits. Plots level curves and percent error in omega.
% This script uses the dielectric function as the "true" solution.
clear; clc;

N = 50; % number of samples in level curves
Nparams = 3;
kappa = 2;

% parameters = [k, sigma, mu]
baseparams = [0.5; 1.5; 900];
txtparams = ["k","\sigma","\mu"];

% set variations
var = 0.50; % x 100% variation considered
xl = (1-var)*baseparams;
xu = (1+var)*baseparams;

% Pre-allocate memory
paramsarr = zeros(Nparams,N);
omega.spectral = zeros(Nparams,N);
omega.xie = zeros(Nparams,N);
omega.xie_rescaled = zeros(Nparams,N);
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

        init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k, sigma, 0, kappa); %tilde{Omega}+igamma

        xi_guess = (init_guess+mu*k)/k;
        xi_guess_rescaled = init_guess/(sigma*k);

        omega.spectral(i,j) = init_guess+mu*k; 
        omega.xie(i,j) = Kappa_Disp_Using_Xie(k,sigma,mu,kappa,xi_guess)*k;
        omega.xie_rescaled(i,j) = Kappa_Disp_Using_Xie(k*sigma,1,0,kappa,xi_guess_rescaled)*sigma*k + mu*k;
        omega.dielectric(i,j) = Kappa_dielectric(k,sigma,0,kappa,init_guess)+mu*k;
    end

    % re-save params as base parameters for next level curve
    params = baseparams;
end

% calculate percent error
errorRe.spectral = abs(real(omega.dielectric)-real(omega.spectral))./abs(real(omega.dielectric));
errorRe.xie = abs(real(omega.dielectric)-real(omega.xie))./abs(real(omega.dielectric));
errorRe.xie_rescaled = abs(real(omega.dielectric)-real(omega.xie_rescaled))./abs(real(omega.dielectric));

errorIm.spectral = abs(imag(omega.dielectric)-imag(omega.spectral))./abs(imag(omega.dielectric));
errorIm.xie = abs(imag(omega.dielectric)-imag(omega.xie))./abs(imag(omega.dielectric));
errorIm.xie_rescaled = abs(imag(omega.dielectric)-imag(omega.xie_rescaled))./abs(imag(omega.dielectric));

%% Plot Level Curves and Percent Errors
close all
txtbase = string(['$(k,\sigma,\mu)$=(',num2str(baseparams(1)),...
    ',',num2str(baseparams(2)),',',num2str(baseparams(3)),')']);

% plot the level curves
figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    plot(paramsarr(i,:), imag(omega.spectral(i,:)),'.-'); hold on
    plot(paramsarr(i,:), imag(omega.xie(i,:)),'.-');
    plot(paramsarr(i,:), imag(omega.xie_rescaled(i,:)),'.-');
    plot(paramsarr(i,:), imag(omega.dielectric(i,:)),'k');
        title('Kappa $\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',16);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',14);
        ylabel('$\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',14);
end
    leg = legend('Spectral','Xie','Rescaled Xie','Dielectric Func.','Orientation','Horizontal');
    leg.Layout.Tile = 'north';
    sgtitle('$\gamma$ Level Curves around '+txtbase,'Interpreter','latex','FontSize',16);

figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    plot(paramsarr(i,:), real(omega.spectral(i,:)),'.-'); hold on
    plot(paramsarr(i,:), real(omega.xie(i,:)),'.-');
    plot(paramsarr(i,:), real(omega.xie_rescaled(i,:)),'.-');
    plot(paramsarr(i,:), real(omega.dielectric(i,:)),'k');
        title('Kappa $\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',16);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',14);
        ylabel('$\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',14);
end
    leg = legend('Spectral','Xie','Rescaled Xie','Dielectric Func.','Orientation','Horizontal');
    leg.Layout.Tile = 'north';
    sgtitle('$\Omega$ Level Curves around '+txtbase,'Interpreter','latex','FontSize',16);

% plot the percent error
figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    semilogy(paramsarr(i,:), errorIm.spectral(i,:),'.-'); hold on
    semilogy(paramsarr(i,:), errorIm.xie(i,:),'.-');
    semilogy(paramsarr(i,:), errorIm.xie_rescaled(i,:),'.-');
        title('\% Error in $\gamma('+txtparams(i)+')$','Interpreter','latex','FontSize',16);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',14);
        ylabel('\% Error','Interpreter','latex','FontSize',14);
        ylim([0,max(1,max(errorIm.xie(i,:)))]); grid on;
end
    leg = legend('Spectral','Xie','Rescaled Xie','Orientation','Horizontal');
    leg.Layout.Tile = 'north';
    sgtitle('\% Error in $\gamma$ around '+txtbase,'Interpreter','latex','FontSize',16);

figure; tiledlayout(2,2);
for i=1:Nparams
    nexttile
    semilogy(paramsarr(i,:), errorRe.spectral(i,:),'.-'); hold on
    semilogy(paramsarr(i,:), errorRe.xie(i,:),'.-');
    semilogy(paramsarr(i,:), errorRe.xie_rescaled(i,:),'.-');
        title('\% Error in $\Omega('+txtparams(i)+')$','Interpreter','latex','FontSize',16);
        xlabel('$'+txtparams(i)+'$','Interpreter','latex','FontSize',14);
        ylabel('\% Error','Interpreter','latex','FontSize',14);
        ylim([0,max(1,max(errorRe.xie(i,:)))]); grid on;
end
    leg = legend('Spectral','Xie','Rescaled Xie','Orientation','Horizontal');
    leg.Layout.Tile = 'north';
    sgtitle('\% Error in $\Omega$ around '+txtbase,'Interpreter','latex','FontSize',16);

pause(0.5)
%% Run simulations N^2 times with randomly selected parameters
rng('shuffle');
parfor j = 1:N^2
    rng(sum(100*clock)+pi*j);
    % Randomly sample parameters within acceptable ranges
    Xs(j,:) = 2*rand(1,Nparams) - 1;
    randparams = 1/2*(diag(xu - xl)*Xs(j,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(randparams(1), randparams(2), 0, kappa);
    xi_guess = init_guess/(randparams(2)*randparams(1));

    spectral = init_guess + randparams(3)*randparams(1);
    xie = Kappa_Disp_Using_Xie(randparams(1)*randparams(2),1,0,kappa,xi_guess)*randparams(1)*randparams(2) + randparams(1)*randparams(3);
    dielectric = Kappa_dielectric(randparams(1),randparams(2),0,kappa,init_guess)+randparams(3)*randparams(1);

    errorRand(j).omega = abs(real(dielectric)-real(xie))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(xie))/abs(imag(dielectric));
    errorRand(j).spectral = abs(real(dielectric)-real(spectral))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(spectral))/abs(imag(dielectric));
end

%% Plot the percent error in the N^2 randomly drawn samples
figure; 
subplot(1,2,1)
semilogy(real([errorRand.omega]),'*');
    title('$\Omega$ \% Error','Interpreter','latex','FontSize',16);
    xlabel('Sample \#','Interpreter','latex','FontSize',14);
    ylabel('$\big|\frac{\Omega_{dielctric}-\Omega_{Xie}}{\Omega_{dielectric}}\big|$','Interpreter','latex','FontSize',14);
    xlim([0,N^2]); grid on
subplot(1,2,2)
semilogy(imag([errorRand.omega]),'*');
    title('$\gamma$ \% Error','Interpreter','latex','FontSize',16);
    xlabel('Sample \#','Interpreter','latex','FontSize',14);
    ylabel('$\big|\frac{\gamma_{dielectric}-\gamma_{Xie}}{\gamma_{dielectric}}\big|$','Interpreter','latex','FontSize',14);
    xlim([0,N^2]); grid on

txt = ['Kappa, ',num2str(var*100),'\% variation on $(k,\sigma,\mu)$=(',...
    num2str(baseparams(1)),',',num2str(baseparams(2)),',',num2str(baseparams(3)),')'];
sgtitle(txt,'Interpreter','latex','FontSize',16);

errorMax_spectral = max(abs([errorRand.spectral]))
errorMax = max(abs([errorRand.omega]))