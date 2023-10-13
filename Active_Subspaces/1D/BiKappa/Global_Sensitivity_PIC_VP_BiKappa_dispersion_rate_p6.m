clear; %close all;
rng('shuffle');
% parpool(4); 
% Initialize algorithm parameters
N = 64;      %Number of samples for each parameter
h = 1e-6;     %Finite difference step size

kappa = 2; % kappa values allowed: kappa > 3/2
beta = 0.9;
mu1 = 818000;

% Pre-allocate memory
growth = zeros(N,1);                        %Output of interest (growth rate)
Nparams = 6;
growth_plus = zeros(N,Nparams);             %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
Xs = zeros(N,Nparams);                      %To save the normalized parameters                

% [k; theta1; theta2; mu1; mu2; beta; kappa];
setvals = [0.5; 1; 1; mu1; mu1+2; beta];

% % plot the distribution
% theta1 = setvals(2);    mu1 = setvals(4);    
% theta2 = setvals(3);    mu2 = setvals(5);
% beta = setvals(6);      
% vplot = linspace(mu1*0.9,mu1*1.1,1000);
% C1=(pi*theta1^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
% C2=(pi*theta2^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
% f0=beta*C1*(1+(vplot-mu1).^2/(kappa-1.5)*theta1^2).^(-kappa) + ...
%    (1-beta)*C2*(1+(vplot-mu2).^2/(kappa-1.5)*theta2^2).^(-kappa);
% figure; plot(vplot,f0)

var = 0.01; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu_1 around 0
if setvals(4) == 0
    xl(4) = -var;
    xu(4) = var;
end

xu(6) = min(xu(6),1); % keep beta <= 1

% Run simulation
tic
parfor jj = 1:N
    rng(sum(100*clock)+pi*jj);
    % Randomly sample parameters within acceptable ranges
    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));

    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    % init_guess = Vlasov_1D_linearized_Steve_v4_Kappa( ...
    %     params(1), params(2), params(3), params(4), params(5), params(6), kappa);
    % % growth(jj) = dielectric_kappa( ...
    % %     params(1), params(2), params(3), params(4), params(5), params(6), kappa-1, init_guess);
    % growth(jj) = Kappa_Bump_Disp_Using_Xie(...
    %     params(1), params(2), params(3), params(4), params(5), params(6), kappa, init_guess); 
    
    % rescaled functions: (gm 10.11.23)
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(params(1),params(2),params(3),0,params(5)-params(4),params(6),kappa); %tilde{Omega}+igamma
    xi_guess = init_guess/(params(2)*params(1));
    omega = Kappa_Bump_Disp_Using_Xie(params(1)*params(2),1,params(3)/params(2),0,(params(5)-params(4))/params(2),params(6),kappa,xi_guess)*params(2)*params(1) + params(4)*params(1); %omega=xi*sigma*k+mu*k
    growth(jj) = imag(omega);
    % error(jj) = abs(omega_exact-growth(jj));
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);

        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        % init_guess_plus = Vlasov_1D_linearized_Steve_v4_Kappa( ...
        %     paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), kappa);
        % % growth_plus(jj, kk) = dielectric_kappa( ...
        % %     paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), kappa-1, init_guess_plus);
        % growth_plus(jj, kk) = Kappa_Bump_Disp_Using_Xie(...
        %     paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), kappa, init_guess_plus);
        
        % rescaled functions: (gm 10.11.23)
        init_guess_plus = Vlasov_1D_linearized_Steve_v4_Kappa(paramsplus(1),paramsplus(2),paramsplus(3),0,paramsplus(5)-paramsplus(4),paramsplus(6),kappa); %tilde{Omega}+igamma
        xi_guess_plus = init_guess_plus/(paramsplus(2)*paramsplus(1));
        omega_plus = Kappa_Bump_Disp_Using_Xie(paramsplus(1)*paramsplus(2),1,paramsplus(3)/paramsplus(2),0,(paramsplus(5)-paramsplus(4))/paramsplus(2),paramsplus(6),kappa,xi_guess_plus)*paramsplus(2)*paramsplus(1) + paramsplus(4)*paramsplus(1); %omega=xi*sigma*k+mu*k
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
save(['Data/kappa',int2str(kappa),'/beta' num2str(beta) '/Dispersion_Rate_Kappa_Bump_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'data_par.mat'])

%exit
delete(gcp('nocreate'))
load train, sound(y,Fs)