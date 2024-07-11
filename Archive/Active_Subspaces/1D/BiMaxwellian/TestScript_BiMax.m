%% Bi-Maxwellian Test Script for Accuracy of Methods
% Runs level curves in each variable and generates N^2 random samples
% inside variation limits. Plots level curves and relative error in omega.
% This script uses the dispersion function as the "true" solution, the
% spectral code for the initial guess, and the Xie code shifted, not scaled.
% can run shifted and scaled though
clear; clc;

N = 100; % number of samples in level curves
Nparams = 6;

% parameters = [k, sigma1, sigma2, mu1, mu2, beta]
baseparams = [0.5; 2; 0.6; 1; 5.5; 0.9];
txtparams = ["k","\sigma_1","\sigma_2","\mu_1","\mu_2","\beta"];

% set variations
var = 0.50; % x 100% variation considered
xl = (1-var)*baseparams;
xu = (1+var)*baseparams;
xu(6) = min(xu(6),0.99); % keep beta < 1
xl(6) = max(xl(6),0.51); % keep beta > 0.5

% Pre-allocate memory
paramsarr = zeros(Nparams,N);
omega.init = zeros(Nparams,N);
omega.xie = zeros(Nparams,N);
omega.xie_shift = zeros(Nparams,N);
omega.xie_shiftscaled = zeros(Nparams,N);
omega.dielectric = zeros(Nparams,N);
errorRand.omega = zeros(1,N^2);
errorRand.init = zeros(1,N^2);

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

        init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma1,sigma2,0,mu2-mu1,beta) + mu1*k;
        % init_guess = BohmGross_BiMax(k,sigma1,sigma2,0,mu2-mu1,beta) + mu1*k;

        xi_guess = init_guess/k;
        xi_guess_shift = init_guess/k-mu1;
        xi_guess_shiftscaled = (init_guess/k-mu1)/sigma1;

        omega.init(i,j) = init_guess; 
        omega.xie(i,j) = BiMaxwellian_Disp_Using_Xie(k,sigma1,sigma2,mu1,mu2,beta,xi_guess)*k;
        omega.xie_shift(i,j) = BiMaxwellian_Disp_Using_Xie(k,sigma1,sigma2,0,mu2-mu1,beta,xi_guess_shift)*k + mu1*k;
        omega.xie_shiftscaled(i,j) = BiMaxwellian_Disp_Using_Xie(k*sigma1,1,sigma2/sigma1,0,(mu2-mu1)/sigma1,beta,xi_guess_shiftscaled)*sigma1*k + mu1*k;
        omega.dielectric(i,j) = BiMax_dielectric([k,sigma1,sigma2,0,mu2-mu1,beta],xi_guess-mu1)*k + mu1*k;
    end

    % re-save params as base parameters for next level curve
    params = baseparams;
end

% calculate relative error
errorRe.init = abs(real(omega.dielectric)-real(omega.init))./abs(real(omega.dielectric));
errorRe.xie = abs(real(omega.dielectric)-real(omega.xie))./abs(real(omega.dielectric));
errorRe.xie_shift = abs(real(omega.dielectric)-real(omega.xie_shift))./abs(real(omega.dielectric));
errorRe.xie_shiftscaled = abs(real(omega.dielectric)-real(omega.xie_shiftscaled))./abs(real(omega.dielectric));

errorIm.init = abs(imag(omega.dielectric)-imag(omega.init))./abs(imag(omega.dielectric));
errorIm.xie = abs(imag(omega.dielectric)-imag(omega.xie))./abs(imag(omega.dielectric));
errorIm.xie_shift = abs(imag(omega.dielectric)-imag(omega.xie_shift))./abs(imag(omega.dielectric));
errorIm.xie_shiftscaled = abs(imag(omega.dielectric)-imag(omega.xie_shiftscaled))./abs(imag(omega.dielectric));

%% Run simulations N^2 times with randomly selected parameters
rng('shuffle');
rng(sum(100*clock));
Xs = 2*rand(N^2,Nparams) - 1; % do sampling in serial

parpool(10);
parfor j = 1:N^2
    % Randomly sample parameters within acceptable ranges
    randparams = 1/2*(diag(xu - xl)*Xs(j,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with randomly drawn parameters
    % init_guess = BohmGross_BiMax(randparams(1),randparams(2),randparams(3),0,randparams(5)-randparams(4),randparams(6));
    init_guess = Vlasov_1D_linearized_Steve_v4(randparams(1),randparams(2),randparams(3),0,randparams(5)-randparams(4),randparams(6));
    xi_guess = init_guess/randparams(1); % shifted (not scaled)
 
    spectral = init_guess + randparams(4)*randparams(1);
    xie = BiMaxwellian_Disp_Using_Xie(randparams(1),randparams(2),randparams(3),0,randparams(5)-randparams(4),randparams(6),xi_guess)*randparams(1) + randparams(4)*randparams(1);
    dielectric = BiMax_dielectric([randparams(1),randparams(2),randparams(3),0,randparams(5)-randparams(4),randparams(6)],xi_guess)*randparams(1) + randparams(4)*randparams(1);
        
    errorRand(j).omega = abs(real(dielectric)-real(xie))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(xie))/abs(imag(dielectric));
    errorRand(j).init = abs(real(dielectric)-real(spectral))/abs(real(dielectric))+1i*abs(imag(dielectric)-imag(spectral))/abs(imag(dielectric));
end
delete(gcp('nocreate'));
errorMax_init = max(abs([errorRand.init]));
errorMax = max(abs([errorRand.omega]));
errorMean = mean(abs([errorRand.omega]));

%% Save the testing data 
save(['Data\TestingDataBiMax_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '_data.mat'])