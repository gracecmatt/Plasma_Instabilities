clear; clc
rng('shuffle');
parpool(12);
% Initialize algorithm parameters
N = 480;            %Number of samples for each parameter
h = 1e-6;           %Finite difference step size

% Pre-allocate memory
growth = zeros(N,1);                     %Output of interest (growth rate)
Nparams = 5;
growth_plus = zeros(N,Nparams);          %Perturbed output of interest
grad_growth = zeros(Nparams,N);          %Gradient of output of interest 
Xs = zeros(N,Nparams);                   %To save the normalized parameters                    

mu = 100;
nu = mu-1;
M = 100;

% vals = [k, sigma, mu, nu, M]
setvals = [1; 1; mu; nu; M];

var = 0.25; % x 100% variation considered
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu around 0
if setvals(3) == 0 
    xl(3) = -var;
    xu(3) = var;
end

% Run simulation
tic
parfor jj = 1:N
    rng(sum(100*clock)+pi*jj);
    % Randomly sample parameters within acceptable ranges
    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    init_guess = Vlasov_1D_linearized_Steve_v4(params(1), params(2), 0, params(4),params(5)) + params(1)*params(3); % + mu*k
    growth(jj) = IncompleteMax_Disp_Using_Xie(params(1), params(2), 0, params(4), params(5), init_guess);
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);                     
        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));

        init_guess_plus = Vlasov_1D_linearized_Steve_v4(paramsplus(1),paramsplus(2),0,paramsplus(4),paramsplus(5)) + paramsplus(1)*paramsplus(3);
        growth_plus(jj, kk) = IncompleteMax_Disp_Using_Xie(paramsplus(1),paramsplus(2),0,paramsplus(4),paramsplus(5),init_guess_plus);
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

%Compute the eigenvalues of C
evalues = diag(S.^2);

% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));

%Save the trial data
save(['Data\Dispersion_Rate_IncompleteMax_nu' num2str(nu) '_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'var_data.mat'])

% exit
delete(gcp('nocreate'))