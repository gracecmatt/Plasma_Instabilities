clear; %close all;
rng('shuffle');
parpool(12); 
% Initialize algorithm parameters
N = 12;      %Number of samples for each parameter
h = 1e-8;     %Finite difference step size

kappa = 1.93; % kappa values allowed: kappa > 3/2
beta = 1; % beta values allowed: beta in [0,1]
sigma =  0.007065164793364;

% Pre-allocate memory
growth = zeros(N,1);                        %Output of interest (growth rate)
Nparams = 7;
growth_plus = zeros(N,Nparams);             %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
Xs = zeros(N,Nparams);                      %To save the normalized parameters 
% aN = zeros(N,1);
% aN_plus = zeros(N,Nparams);

% [k; theta_1; theta_2; mu_1; mu_2; beta; kappa];
setvals = [0.9; sigma; 1; 818; 1; beta; kappa];

% % plot the distribution
% theta1 = setvals(2);    mu1 = setvals(4);  
% theta2 = setvals(3);    mu2 = setvals(5);
% beta = setvals(6);
% vplot = -2:0.01:4;
% C1=(pi*theta1^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
% C2=(pi*theta2^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
% f0=beta*C1*(1+(vplot-mu1).^2/((kappa-1.5)*theta1^2)).^(-kappa) + ...
%    (1-beta)*C2*(1+(vplot-mu2).^2/((kappa-1.5)*theta2^2)).^(-kappa);
% plot(vplot,f0,'LineWidth',2)
% ylim([-0.1,2.6])
% txt = sprintf('Fast Bulk Kappa Equilibrium \\kappa=1.9, \\sigma=0.3487');
% title(txt)

var = 0.0; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu_1 around 0
xl(4) = xl(5)-setvals(5);
xu(4) = xu(5)-setvals(5);

xu(6) = min(xu(6),1); % keep beta <= 1
xl(7) = max(xl(7),1.5 + 1e-10); % keep kappa > 3/2

% Run simulation
tic
for jj = 1:N
    rng(sum(100*clock)+pi*jj);
    % Randomly sample parameters within acceptable ranges
    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa( ...
        params(1), params(2), params(3), params(4), params(5), params(6), params(7));
    % growth(jj) = dielectric_kappa( ...
    %     params(1), params(2), params(3), params(4), params(5), params(6), kappa-1, init_guess);
    growth(jj) = Kappa_Bump_Disp_Using_Xie(...
        params(1), params(2), params(3), params(4), params(5), params(6), params(7), init_guess); 
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);

        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        init_guess_plus = Vlasov_1D_linearized_Steve_v4_Kappa( ...
            paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), paramsplus(7));
        % growth_plus(jj, kk) = dielectric_kappa( ...
        %     paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), kappa-1, init_guess_plus);
        growth_plus(jj,kk) = Kappa_Bump_Disp_Using_Xie(...
            paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), paramsplus(7), init_guess_plus);
    end
end
parfor jj = 1:N
    % Calculate the appx gradients using finite differences
    grad_growth(:,jj) = (growth_plus(jj, :) - growth(jj))/h;
end
toc

% Compute the singular value decomposition of the matrix C
[U,S,V] = svd(1/sqrt(N)*grad_growth);
w = U(:,1);
w2 = U(:,2);

% Compute the eigenvalues of C
evalues = diag(S.^2);

% Compute the condition number (added by gm 8/9/23)
cond = evalues(1)/sum(evalues);
% aN_mean = mean([aN,aN_plus],'all');
% aN_max = max([aN,aN_plus],[],'all');

% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));

%Save the trial data
% save(['Data/kappa',int2str(kappa),'/beta' num2str(beta) '/Dispersion_Rate_Kappa_Bump_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'data_par_2.mat'])
save(['Data/Maksimovi1997/beta' num2str(beta) '/Dispersion_Rate_Kappa_Bump_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'data_par_SLOW.mat'])

%exit
delete(gcp('nocreate'))