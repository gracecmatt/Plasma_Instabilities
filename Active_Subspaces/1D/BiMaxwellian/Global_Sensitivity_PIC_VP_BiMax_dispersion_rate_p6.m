clear; clc;
rng('shuffle');
parpool(8);
% Initialize algorithm parameters
N = 512;                              %Number of samples for each parameter
h = 1e-6;                                      %Finite difference step size

% Pre-allocate memory
growth = zeros(N,1);                      %Output of interest (growth rate)
Nparams = 6;
growth_plus = zeros(N,Nparams);               %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
Xs = zeros(N,Nparams);                   %To save the normalized parameters
omega_error = zeros(N,1);

% vals = [k, sigma1, sigma2, mu1, mu2, beta]
setvals = [0.5; 2; 0.6; 1; 5.5; 0.9];

var = 0.50; % x 100% variation considered
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu1 around 0
if setvals(4)==0
    xl(4) = -var;
    xu(4) = var;
end

xu(6) = min(xu(6),0.97); % keep beta < 1
xl(6) = max(xl(6),0.58); % keep beta > 0.5

% Run simulation
tic
parfor jj = 1:N
    rng(sum(100*clock)+pi*jj);
    % Randomly sample parameters within acceptable ranges
    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with randomly drawn parameters
    init_guess = Vlasov_1D_linearized_Steve_v4(params(1),params(2),params(3),0,params(5)-params(4),params(6));
    xi_guess = init_guess/params(1); % shifted (not scaled)
 
    omega = BiMaxwellian_Disp_Using_Xie(params(1),params(2),params(3),0,params(5)-params(4),params(6),xi_guess)*params(1) + params(4)*params(1);
    dielectric = BiMax_dielectric([params(1),params(2),params(3),0,params(5)-params(4),randparams(6)],xi_guess)*params(1) + params(4)*params(1);
    growth(jj) = imag(omega);
    omega_error(jj) = abs(real(omega)-real(dielectric)) + 1i*abs(imag(omega)-imag(dielectric));
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);                     
        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));

        init_guess_plus = Vlasov_1D_linearized_Steve_v4(paramsplus(1),paramsplus(2),paramsplus(3),0,paramsplus(5)-paramsplus(4),paramsplus(6));
        xi_guess_plus = init_guess_plus/paramsplus(1); % shifted (not scaled)
     
        omega0_plus = BiMaxwellian_Disp_Using_Xie(paramsplus(1),paramsplus(2),paramsplus(3),0,paramsplus(5)-paramsplus(4),paramsplus(6),xi_guess_plus)*paramsplus(1);
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

%Compute the eigenvalues of C
evalues = diag(S.^2);

% Compute the condition number
cond = evalues(1)/sum(evalues);

% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));

%Save the trial data
save(['Data\Dispersion_BiMax_P' int2str(Nparams) '_N' int2str(N) '_var' num2str(var) '_data.mat'])

% exit
delete(gcp('nocreate'))