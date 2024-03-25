clear; clc;
rng('shuffle');
parpool(12); 
% Initialize algorithm parameters
N = 504;                              %Number of samples for each parameter
h = 1e-6;                                      %Finite difference step size

kappa = 4;                              % kappa values allowed: kappa > 3/2              

% Pre-allocate memory
Nparams = 6;
growth = zeros(N,1);                      %Output of interest (growth rate)
growth_plus = zeros(N,Nparams);               %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest               
omega_error = zeros(N,1);           %Compare Xie/dielectric funtion results

% vals = [k; sigma1; sigma2; mu; v0; beta];
setvals = [0.5; 1; 0.6; 1; 5; 0.95];

var = 0.15; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu_1 around 0
if setvals(4) == 0
    xl(4) = -var;
    xu(4) = var;
end

xu(6) = min(xu(6),0.99); % keep beta < 1
xl(6) = max(xl(6),0.55); % keep beta > 0.5

% Run simulation
tic
rng(sum(100*clock));
Xs = 2*rand(N,Nparams) - 1; % do sampling in serial
parfor jj = 1:N
    % Randomly sample parameters within acceptable ranges
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));

    % rescaled functions: (gm 10.11.23)
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(params(1),params(2),params(3),0,params(5),params(6),kappa);
    xi_guess = init_guess/(params(2)*params(1)); % shifted and scaled

    omega = BiKappa_Disp_Using_Xie(params(1)*params(2),1,params(3)/params(2),0,params(5)/params(2),params(6),kappa,xi_guess)*params(2)*params(1) + params(4)*params(1);
    dielectric = BiKappa_dielectric(params(1),params(2),params(3),0,params(5),params(6),kappa,init_guess) + params(4)*params(1);
    growth(jj) = imag(omega);
    omega_error(jj) = abs(real(omega)-real(dielectric)) + 1i*abs(imag(omega)-imag(dielectric));
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));

        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        init_guess_plus = Vlasov_1D_linearized_Steve_v4_Kappa(paramsplus(1),paramsplus(2),paramsplus(3),0,paramsplus(5),paramsplus(6),kappa);
        xi_guess_plus = init_guess_plus/(paramsplus(2)*paramsplus(1));
        
        omega0_plus = BiKappa_Disp_Using_Xie(paramsplus(1)*paramsplus(2),1,paramsplus(3)/paramsplus(2),0,paramsplus(5)/paramsplus(2),paramsplus(6),kappa,xi_guess_plus)*paramsplus(2)*paramsplus(1);
        growth_plus(jj, kk) = imag(omega0_plus);
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

% Compute the condition number
cond = evalues(1)/sum(evalues);

% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));

%Save the trial data
save(['Dispersion_BiKappa' int2str(kappa) '_P' int2str(Nparams) '_N' int2str(N) '_var' int2str(var*100) '_data.mat'])

%exit
delete(gcp('nocreate'))