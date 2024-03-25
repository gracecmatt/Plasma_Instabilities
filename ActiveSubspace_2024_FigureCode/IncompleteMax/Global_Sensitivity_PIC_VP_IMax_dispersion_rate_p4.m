clear; clc
rng('shuffle');
parpool(12);
% Initialize algorithm parameters
N = 1200;            %Number of samples for each parameter
h = 1e-6;           %Finite difference step size

% Pre-allocate memory
Nparams = 4;
growth = zeros(N,1);                      %Output of interest (growth rate)
growth_plus = zeros(N,Nparams);               %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest                 

% vals = [k, sigma, mu, nu]
setvals = [0.5; 1; 0; -1];

var = 0.5; % x 100% variation considered
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu around 0
if setvals(3) == 0 
    xl(3) = -var;
    xu(3) = var;
end

% Run simulation
tic
rng(sum(100*clock));
Xs = 2*rand(N,Nparams) - 1; % do sampling in serial
parfor jj = 1:N
    % Randomly sample parameters within acceptable ranges
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    init_guess = BohmGross_IMax(params(1),params(2),0,params(4));
    xi_guess = init_guess/(params(2)*params(1)); % shifted & scaled
    
    omega = IncompleteMax_Disp_Using_Xie(params(1)*params(2),1,0,(params(4)-params(3))/params(2),xi_guess)*params(2)*params(1);
    growth(jj) = imag(omega);
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);                     
        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));

        init_guess_plus = BohmGross_IMax(paramsplus(1),paramsplus(2),0,paramsplus(4));
        xi_guess_plus = init_guess_plus/(paramsplus(2)*paramsplus(1)); % shifted & scaled

        omega_plus = IncompleteMax_Disp_Using_Xie(paramsplus(1)*paramsplus(2),1,0,(paramsplus(4)-paramsplus(3))/paramsplus(2),xi_guess_plus)*paramsplus(2)*paramsplus(1);
        growth_plus(jj,kk) = imag(omega_plus);
    end
end

parfor jj = 1:N
    % Calculate the appx gradients using finite differences
    grad_growth(:,jj) = (growth_plus(jj, :) - growth(jj))/h;
end
toc; sound(sin(1:.4:400))

% Compute the singular value decomposition of the matrix C
[U,S,V] = svd(1/sqrt(N)*grad_growth);
w = U(:,1);
w2 = U(:,2);

%Compute the eigenvalues of C
evalues = diag(S.^2);

% Compute the first two eigenvalue ratios
eta(1) = evalues(1)/sum(evalues);
eta(2) = (evalues(1)+evalues(2))/sum(evalues);

% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));

%Save the trial data
save(['Data/Dispersion_IncompleteMax_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var*100) '_data.mat'])

% exit
delete(gcp('nocreate'))