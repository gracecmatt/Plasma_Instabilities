clear; %close all;
rng('shuffle');
parpool(8); 
% Initialize algorithm parameters
N = 512*4;      %Number of samples for each parameter
h = 1e-6;     %Finite difference step size

kappa = 2; % kappa values allowed: kappa > 3/2
beta = 0.9;
mu1 = 0;

% Pre-allocate memory
growth = zeros(N,1);                        %Output of interest (growth rate)
Nparams = 6;
growth_plus = zeros(N,Nparams);             %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
Xs = zeros(N,Nparams);                      %To save the normalized parameters                
omega_error = zeros(N,1);
spectral_error = zeros(N,1);

% [k; sigma1; sigma2; mu1; mu2; beta; kappa];
setvals = [0.5; 1; 0.6; mu1; mu1+2.5; beta];

var = 0.5; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu_1 around 0
if setvals(4) == 0
    xl(4) = -var;
    xu(4) = var;
end

xu(6) = min(xu(6),0.98); % keep beta <= 1

% Run simulation
tic
parfor jj = 1:N
    rng(sum(100*clock)+pi*jj);
    % Randomly sample parameters within acceptable ranges
    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));

    % rescaled functions: (gm 10.11.23)
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(params(1),params(2),params(3),0,params(5)-params(4),params(6),kappa); %tilde{Omega}+igamma
    omega_guess = init_guess+params(4)*params(1);
    xi_guess = init_guess/(params(2)*params(1));
    omega = BiKappa_Disp_Using_Xie(params(1)*params(2),1,params(3)/params(2),0,(params(5)-params(4))/params(2),params(6),kappa,xi_guess)*params(2)*params(1) + params(4)*params(1); %omega=xi*sigma*k+mu*k
    omega_disp = BiKappa_dielectric(params(1),params(2),params(3),params(4),params(5),params(6),kappa,omega_guess);
    growth(jj) = imag(omega);
    omega_error(jj) = abs(real(omega_disp)-real(omega))+1i*abs(imag(omega_disp)-imag(omega));
    spectral_error(jj) = abs(real(omega_disp)-real(init_guess+params(4)*params(1)))+1i*abs(imag(omega_disp)-imag(init_guess));
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));

        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        init_guess_plus = Vlasov_1D_linearized_Steve_v4_Kappa(paramsplus(1),paramsplus(2),paramsplus(3),0,paramsplus(5)-paramsplus(4),paramsplus(6),kappa); %tilde{Omega}+igamma
        xi_guess_plus = init_guess_plus/(paramsplus(2)*paramsplus(1));
        omega_plus = BiKappa_Disp_Using_Xie(paramsplus(1)*paramsplus(2),1,paramsplus(3)/paramsplus(2),0,(paramsplus(5)-paramsplus(4))/paramsplus(2),paramsplus(6),kappa,xi_guess_plus)*paramsplus(2)*paramsplus(1) + paramsplus(4)*paramsplus(1); %omega=xi*sigma*k+mu*k
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

% Compute the eigenvalues of C
evalues = diag(S.^2);

% Compute the condition number
cond = evalues(1)/sum(evalues);

% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));

%Save the trial data
save(['Data/kappa',int2str(kappa),'/Dispersion_Rate_Kappa_Bump_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'data_par.mat'])

%exit
delete(gcp('nocreate'))
% sound(sin(1:3000));

% plot the error
figure; 
subplot(1,2,1)
semilogy(real(omega_error),'*'); grid on
% ylim([10^(-17),10^(-5)]);
xlabel('Sample \#','Interpreter','latex','FontSize',12)
ylabel('$|\Omega_{dielectirc}-\Omega_{Xie}|$','Interpreter','latex','FontSize',12)
title('$\Omega$ Error','Interpreter','latex','FontSize',12);
subplot(1,2,2)
semilogy(imag(omega_error),'*'); grid on
% ylim([10^(-17),10^(-5)]);
xlabel('Sample \#','Interpreter','latex','FontSize',12);
ylabel('$|\gamma_{dielectric}-\gamma_{Xie}|$','Interpreter','latex','FontSize',12)
title('$\gamma$ Error','Interpreter','latex','FontSize',12);
txt = [num2str(var*100),'\% variation on $(k,\sigma_1,\sigma_2,\mu_1,\mu_2,\beta)$=(',num2str(setvals(1)),',',num2str(setvals(2)),',',num2str(setvals(3)),',',num2str(setvals(4)),',',num2str(setvals(5)),',',num2str(setvals(6)),')'];
sgtitle(txt,'Interpreter','latex','FontSize',14)

max_error = max(abs(omega_error))

max_spectral_error = max(abs(spectral_error))