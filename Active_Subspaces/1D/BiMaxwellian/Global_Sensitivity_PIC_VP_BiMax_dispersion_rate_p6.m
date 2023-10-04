clear; clc;
rng('shuffle');
parpool(6); %was 32
% Initialize algorithm parameters
N = 240;                              %Number of samples for each parameter
h = 1e-6;                                      %Finite difference step size

% Pre-allocate memory
growth = zeros(N,1);                      %Output of interest (growth rate)
Nparams = 6;
growth_plus = zeros(N,Nparams);               %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
Xs = zeros(N,Nparams);                   %To save the normalized parameters
%w = zeros(Nparams,1);                                %Weight vectors 
%evalues = zeros(Nparams,1);           %Eigenvalues of the C matrix
%diff_growth = 0; %zeros(1,1);              %Differences in largest and
%smallest element of grad_growth
%I = eye(Nparams);                     

% vals = [k, sigma1, sigma2, mu1, mu2, beta]
setvals = [0.5; 1; 1; 0; 2.5; 0.9];

var = 0.15; % x 100% variation considered
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu1 around 0
if setvals(4)==0
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
    k = params(1);
    sigma1 = params(2);
    sigma2 = params(3);
    mu1 = params(4);
    mu2 = params(5);
    beta = params(6);

    % % % Numerically solve 1D Vlasov-Poisson with baseline parameters
    % % init_guess = Vlasov_1D_linearized_Steve_v4(params(1),params(2),params(3),0,params(5)-params(4),params(6)); %tilde{Omega}+igamma
    % % xi_guess = init_guess/(params(1)*params(2));
    % % omega = BiMaxwellian_Disp_Using_Xie(params(1)*params(2),1,params(3)/params(2),0,...
    % %     (params(5)-params(4))/params(2),params(6),xi_guess)*params(1)*params(2) + params(1)*params(4); %omega=xi*k*sigma1+mu1*k
    % % growth(jj) = imag(omega);

    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma1,sigma2,0,mu2-mu1,beta); %tilde{Omega}+igamma
    xi_guess = (init_guess)/(k);
    % omega = BiMaxwellian_Disp_Using_Xie(k, sigma1, sigma2, mu1, mu2, beta, xi_guess)*k; %omega=xi*k
    omega = BiMaxwellian_Disp_Using_Xie(k*sigma1,1,sigma2/sigma1,0,(mu2-mu1)/sigma1,beta,xi_guess)*k*sigma1 + k*mu1; %omega=xi*k*sigma1+mu1*k
    growthZetaf =  imag(dispersion_growthrate_BiMax(params,(init_guess+k*mu1)/k))*k;
    growth(jj) = imag(omega);
    error(jj) = abs(growthZetaf-imag(omega));

    % % init_guess = Vlasov_1D_linearized_Steve_v4(params(1), params(2), 0);
    % % xi_guess = init_guess/(params(1)*params(2));
    % % omega(jj) = Maxwellian_Disp_Using_Xie(params(1)*params(2), 1, 0, xi_guess)*params(1)*params(2) + params(1)*params(3);
    % % growth(jj) = imag(omega(jj));

    % growth(jj) = dispersion_growthrate_BiMax(params, init_guess);
    
    %while growth(jj) < 0 || growth(jj) > 2 
    % while (growth(jj) > 5  || growth(jj) < 1e-10)
    %     Xs(jj,:) = 2*rand(1,Nparams) - 1;
    %     params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
    %     growth(jj) = dispersion_growthrate_BiMax(params);
    % end 

    % % rescaled version:
    % init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma1, sigma2, 0, mu2-mu1, beta); %\tilde{Omega}+igamma
    % xi_scaled = init_guess/(sigma1*k);
    % omega_xie_rescaled(count) = BiMaxwellian_Disp_Using_Xie(k*sigma1, 1, sigma2/sigma1, 0, (mu2-mu1)/sigma1, beta, xi_scaled)*k*sigma1 + mu1*k; % omega = xi*k*sigma1 + mu1*k
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);                     
        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));

        k = paramsplus(1);
        sigma1 = paramsplus(2);
        sigma2 = paramsplus(3);
        mu1 = paramsplus(4);
        mu2 = paramsplus(5);
        beta = paramsplus(6);

        % % init_guess_plus = Vlasov_1D_linearized_Steve_v4(paramsplus(1),paramsplus(2),paramsplus(3),0,paramsplus(5)-paramsplus(4),paramsplus(6));
        % % xi_guess_plus = init_guess_plus/(paramsplus(1)*paramsplus(2));
        % % omega_plus = BiMaxwellian_Disp_Using_Xie(paramsplus(1)*paramsplus(2),1,paramsplus(3)/paramsplus(2),0,...
        % %     (paramsplus(5)-paramsplus(4))/paramsplus(2),paramsplus(6),xi_guess_plus)*paramsplus(1)*paramsplus(2) + paramsplus(1)*paramsplus(4);
        % % growth_plus(jj, kk) = imag(omega_plus);
        % % % growth_plus(jj, kk) = dispersion_growthrate_BiMax(paramsplus,init_guess_plus);

        % Numerically solve 1D Vlasov-Poisson with baseline parameters
        init_guess_plus = Vlasov_1D_linearized_Steve_v4(k,sigma1,sigma2,0,mu2-mu1,beta); %tilde{Omega}+igamma
        xi_guess_plus = (init_guess_plus)/(k);
        % omega_plus = BiMaxwellian_Disp_Using_Xie(k, sigma1, sigma2, mu1, mu2, beta, xi_guess_plus)*k; %omega=xi*k
        omega_plus = BiMaxwellian_Disp_Using_Xie(k*sigma1,1,sigma2/sigma1,0,(mu2-mu1)/sigma1,beta,xi_guess_plus)*k*sigma1 + k*mu1; %omega=xi*k*sigma1+mu1*k
        % growthZetaf =  imag(dispersion_growthrate_BiMax(paramsplus,(init_guess_plus+mu1*k)/k))*k;
        growth_plus(jj, kk) = imag(omega_plus);
        
        %while growth_plus(jj,kk) <0 || growth_plus(jj,kk) >2 % set to 2 for 10%, set to 1 for 25% runs
        % while (growth_plus(jj,kk) > 5 || growth_plus(jj,kk) < 1e-10)
        %     xplus = randparams + h*I(:,kk);
        %     paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        %     growth_plus(jj, kk)= dispersion_growthrate_BiMax(paramsplus);
        % end 
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
% save(['Data\Dispersion_Rate_BiMax_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'data_par.mat'])

% exit
delete(gcp('nocreate'))