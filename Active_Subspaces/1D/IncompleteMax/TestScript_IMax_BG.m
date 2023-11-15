%% Incomplete Maxwellian Test Script for Accuracy of Methods
% Runs level curves in each variable and generates N^2 random samples
% inside variation limits. Plots level curves in omega.
% This script uses the Maxwellian as the "true" solution for small enough nu, 
% the spectral code for the initial guess, and the Xie code shifted and scaled.
clear; clc;

N = 40; % number of samples in level curves
Nparams = 4;
M = 9; % slope of erf approximation to step function, for spectral code

% parameters = [k, sigma, mu, nu]
baseparams = [0.5; 2; 1; -20];
txtparams = ["k","\sigma","\mu","\nu"];

% set variations
var = 0.50; % x 100% variation considered,
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
            nu=params(4);
    
        init_guess = BohmGross_IMax(k,sigma,0,nu) + mu*k;
%         init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0, nu, M) + mu*k;

        xi_guess = init_guess/k;
        xi_guess_shift = init_guess/k-mu;
        xi_guess_shiftscaled = (init_guess/k-mu)/sigma;

        omega.spectral(i,j) = init_guess; 
        omega.xie(i,j) = IncompleteMax_Disp_Using_Xie(k,sigma,mu,nu,xi_guess)*k;
        omega.xie_shift(i,j) = IncompleteMax_Disp_Using_Xie(k,sigma,0,(nu-mu),xi_guess_shift)*k + mu*k;
        omega.xie_shiftscaled(i,j) = IncompleteMax_Disp_Using_Xie(sigma*k,1,0,(nu-mu)/sigma,xi_guess_shiftscaled)*k*sigma + mu*k;
        omega.dielectric(i,j) = Max_dielectric([k,sigma,0],xi_guess-mu)*k + mu*k;
    end

    % re-save params as base parameters for next level curve
    params = baseparams;
end

% % calculate percent error
errorRe.spectral = abs(real(omega.dielectric)-real(omega.spectral))./abs(real(omega.dielectric));
errorRe.xie = abs(real(omega.dielectric)-real(omega.xie))./abs(real(omega.dielectric));
errorRe.xie_shift = abs(real(omega.dielectric)-real(omega.xie_shift))./abs(real(omega.dielectric));
errorRe.xie_shiftscaled = abs(real(omega.dielectric)-real(omega.xie_shiftscaled))./abs(real(omega.dielectric));

errorIm.spectral = abs(imag(omega.dielectric)-imag(omega.spectral))./abs(imag(omega.dielectric));
errorIm.xie = abs(imag(omega.dielectric)-imag(omega.xie))./abs(imag(omega.dielectric));
errorIm.xie_shift = abs(imag(omega.dielectric)-imag(omega.xie_shift))./abs(imag(omega.dielectric));
errorIm.xie_shiftscaled = abs(imag(omega.dielectric)-imag(omega.xie_shiftscaled))./abs(imag(omega.dielectric));

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
    init_guess = BohmGross_IMax(randparams(1),randparams(2),0,randparams(4));
    % init_guess = Vlasov_1D_linearized_Steve_v4(randparams(1),randparams(2),0,randparams(4),M);
    xi_guess = init_guess/(randparams(2)*randparams(1));
 
    spectral = init_guess + randparams(3)*randparams(1);
    xie = IncompleteMax_Disp_Using_Xie(randparams(1)*randparams(2),1,0,(randparams(4)-randparams(3))/randparams(2),xi_guess)*randparams(2)*randparams(1) + randparams(3)*randparams(1);
    dielectric = Max_dielectric([randparams(1),randparams(2),0],init_guess/randparams(1))*randparams(1) + randparams(3)*randparams(1);
        
    errorRand(j).omega = abs(real(dielectric)-real(xie))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(xie))/abs(imag(dielectric));
    errorRand(j).spectral = abs(real(dielectric)-real(spectral))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(spectral))/abs(imag(dielectric));
end
delete(gcp('nocreate'));
errorMax_spectral = max(abs([errorRand.spectral]));
errorMax = max(abs([errorRand.omega]));
errorMean = mean(abs([errorRand.omega]));
% toc

%% Save the testing data and figures
% eps is a vector format file type that works well in LaTeX
% pdf is a vector format file type that is easy to share and view

save(['TestingDataIMax_BG_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '_data.mat'])

% savefig(fig1,  ['Figs\LevelCurvesIMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig1, ['Figs\LevelCurvesIMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig1, ['Figs\LevelCurvesIMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% savefig(fig2,  ['Figs\LevelCurvesIMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig2, ['Figs\LevelCurvesIMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig2, ['Figs\LevelCurvesIMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% savefig(fig3,  ['Figs\RelErrorIMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig3, ['Figs\RelErrorIMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig3, ['Figs\RelErrorIMax_gamma_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% savefig(fig4,  ['Figs\RelErrorIMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.fig'])
% hgexport(fig4, ['Figs\RelErrorIMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig4, ['Figs\RelErrorIMax_Omega_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% savefig(fig5,  ['Figs\SampRandIMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.fig'])
% hgexport(fig5, ['Figs\SampRandIMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
% hgexport(fig5, ['Figs\SampRandIMax_P' int2str(Nparams) '_N' int2str(N^2) '_var' num2str(var*100) '.pdf'], hgexport('factorystyle'), 'Format', 'pdf');
% 
% clear fig1 fig2 fig3 fig4 fig5 % comment out to save data for figures too
% save(['Data\TestingDataIMax_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '_data.mat'])