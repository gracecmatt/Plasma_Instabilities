rng('shuffle');
% parpool(4);
parpool(32); %was 32
% Initialize algorithm parameters
N = 512;     %Number of samples for each parameter
h = 1e-6;    %Finite difference step size
%trial = 1;   %Trial number (used when saving figures)

kappa = 1;

% Pre-allocate memory
growth = zeros(N,1);                        %Output of interest (growth rate)
Nparams = 6;
growth_plus = zeros(N,Nparams);             %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
Xs = zeros(N,Nparams);                    %To save the normalized parameters
%w = zeros(Nparams,1);                 %Weight vectors 
%evalues = zeros(Nparams,1);           %Eigenvalues of the C matrix
%diff_growth = 0; %zeros(1,1);              %Differences in largest and
%smallest element of grad_growth
%I = eye(Nparams);                     

% [k; theta_1; theta_2; mu_1; mu_2; beta]; 
setvals = [0.5; 1; 1; 0; 4; 0.8];

var = 0.01; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fixing mu_1
xl(4) = xl(5)-4;
xu(4) = xu(5)-4;

% Run simulation
tic
parfor jj = 1:N
    rng(sum(100*clock)+pi*jj);
    % Randomly sample parameters within acceptable ranges
    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa( ...
        params(1), params(2), params(3), params(4), params(5), params(6), kappa);
    growth(jj) = dielectric_kappa( ... % this function assumes kappa=1
        params(1), params(2), params(3), params(4), params(5), params(6), kappa, init_guess);
    %while growth(jj) < 0 || growth(jj) > 2 
    % Will get trapped in this while loop if gamma is always negative
    %while (growth(jj) > 5  || growth(jj) < 1e-10)
    %    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    %    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
    %    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa( ...
    %        params(1), params(2), params(3), params(4), params(5), params(6), kappa);
    %    growth(jj) = Kappa_Bump_Disp_Using_Xie(...
    %        params(1), params(2), params(3), params(4), params(5), params(6), kappa, init_guess);
    %end 
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);                     
        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        init_guess_plus = Vlasov_1D_linearized_Steve_v4_Kappa( ...
            paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), kappa);
        growth_plus(jj) = dielectric_kappa( ...
            paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), kappa, init_guess_plus);
        %while growth_plus(jj,kk) <0 || growth_plus(jj,kk) >2 % set to 2 for 10%, set to 1 for 25% runs
        %while (growth_plus(jj,kk) > 5 || growth_plus(jj,kk) < 1e-10)
        %    xplus = randparams + h*I(:,kk);
        %    paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        %    init_guess_plus = Vlasov_1D_linearized_Steve_v4_Kappa(...
        %        paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), kappa);
        %    growth_plus(jj, kk) = Kappa_Bump_Disp_Using_Xie(...
        %        paramsplus(1), paramsplus(2), paramsplus(3), paramsplus(4), paramsplus(5), paramsplus(6), kappa, init_guess);
        %end 
    end
end
parfor jj = 1:N
    % Calculate the appx gradients using finite differences
    grad_growth(:,jj) = (growth_plus(jj, :) - growth(jj))/h;
end
toc
% Compute the weights, eigenvalues, and plot results

% Compute the singular value decomposition of the matrix C
[U,S,V] = svd(1/sqrt(N)*grad_growth);
w = U(:,1);
w2 = U(:,2);
    
%Compute the eigenvalues of C
evalues = diag(S.^2);
    
% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));
   
%Save the trial data
save(['Dispersion_Rate_Kappa_Bump_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'data_par.mat'])

%exit
