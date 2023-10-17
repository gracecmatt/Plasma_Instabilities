% clear; clc
rng('shuffle');
parpool(4); %was 32
% Initialize algorithm parameters
N = 256;            %Number of samples for each parameter
h = 1e-6;           %Finite difference step size
%trial = 1;         %Trial number (used when saving figures)

% Pre-allocate memory
growth = zeros(N,1);                     %Output of interest (growth rate)
Nparams = 3;
growth_plus = zeros(N,Nparams);          %Perturbed output of interest
grad_growth = zeros(Nparams,N);          %Gradient of output of interest 
Xs = zeros(N,Nparams);                   %To save the normalized parameters
%w = zeros(Nparams,1);                   %Weight vectors 
%evalues = zeros(Nparams,1);           %Eigenvalues of the C matrix
%diff_growth = 0; %zeros(1,1);              %Differences in largest and
%smallest element of grad_growth
%I = eye(Nparams);                     
omega_error = zeros(N,1);   

% vals = [k, sigma, mu]
setvals = [0.5; 5; 100];

var = 0.10; % x 100% variation considered
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
    init_guess = Vlasov_1D_linearized_Steve_v4(params(1), params(2), 0);
    xi_guess = init_guess/(params(2)*params(1));
    omega = Lorentzian_Disp_Using_Xie(params(1)*params(2), 1, 0, xi_guess)*params(1)*params(2) + params(1)*params(3);
    omega_exact = params(3)*params(1)+1 + 1i*(-params(2)*params(1));

    growth(jj) = imag(omega);
    omega_error(jj) = abs(real(omega_exact)-real(omega))+1i*abs(imag(omega_exact)-imag(omega));
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        
        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        init_guess_plus = Vlasov_1D_linearized_Steve_v4(paramsplus(1),paramsplus(2),0);
        xi_guess_plus = init_guess_plus/(paramsplus(2)*paramsplus(1))
        omega_plus = Lorentzian_Disp_Using_Xie(paramsplus(1)*paramsplus(2),1,0,xi_guess_plus)*paramsplus(1)*paramsplus(2) + paramsplus(1)*paramsplus(3);
        growth_plus(jj, kk) = imag(omega_plus);
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

% Compute the condition number (added by gm 10/5/23)
cond = evalues(1)/sum(evalues);

% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));

%Save the trial data
save(['Data\Dispersion_Rate_Lorentz_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'var_data.mat'])

% exit
delete(gcp('nocreate'))

% plot the error
figure; 
subplot(1,2,1)
semilogy(real(omega_error),'*'); grid on
title('$\Omega$ Error','Interpreter','latex','FontSize',12);
subplot(1,2,2)
semilogy(imag(omega_error),'*'); grid on
title('$\gamma$ Error','Interpreter','latex','FontSize',12);
max_error = max(abs(omega_error))