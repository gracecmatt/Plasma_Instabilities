rng('shuffle');
parpool(32);%parpool(24);%was 32
% Initialize algorithm parameters
N = 512;     %Number of samples for each parameter
h = 1e-6;    %Finite difference step size
%trial = 1;   %Trial number (used when saving figures)

kappa = 1;

% Pre-allocate memory
growth = zeros(N,1);                        %Output of interest (growth rate)
Nparams = 3;
growth_plus = zeros(N,Nparams);             %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
Xs = zeros(N,Nparams);                    %To save the normalized parameters
%w = zeros(Nparams,1);                 %Weight vectors 
%evalues = zeros(Nparams,1);           %Eigenvalues of the C matrix
%diff_growth = 0; %zeros(1,1);              %Differences in largest and
%smallest element of grad_growth
%I = eye(Nparams);                     

% Set upper and lower bounds for parameters - beta_2 = 1 - beta_1
% [k, theta, mu]
setvals = [0.5; 1; 0];

var = 0.50; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fixing mu
xl(3) = -2;
xu(3) = 2;

% Run simulation
tic
parfor jj = 1:N
    rng(sum(100*clock)+pi*jj);
    % Randomly sample parameters within acceptable ranges
    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    init_guess = Vlasov_1D_linearized_Steve_v3_Kappa(params(1), params(2), params(3), kappa);
    growth(jj) = Kappa_Disp_Using_Xie(params(1), params(2), params(3), kappa, init_guess);
    %while growth(jj) < 0 || growth(jj) > 2 
    % Will get trapped in this while loop if gamma is always negative
    %while (growth(jj) > 5  || growth(jj) < 1e-10)
    %    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    %    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
    %    init_guess = Vlasov_1D_linearized_Steve_v3_Kappa(params(1), params(2), kappa);
    %    growth(jj) = Kappa_Disp_Relation(params(1), params(2), init_guess);
    %end 
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);                     
        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        init_guess_plus = Vlasov_1D_linearized_Steve_v3_Kappa(paramsplus(1), paramsplus(2), paramsplus(3), kappa);
        growth_plus(jj, kk) = Kappa_Disp_Using_Xie(paramsplus(1), paramsplus(2), paramsplus(3), kappa, init_guess_plus);
        %while growth_plus(jj,kk) <0 || growth_plus(jj,kk) >2 % set to 2 for 10%, set to 1 for 25% runs
        %while (growth_plus(jj,kk) > 5 || growth_plus(jj,kk) < 1e-10)
        %   xplus = randparams + h*I(:,kk);
        %   paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        %   init_guess_plus = Vlasov_1D_linearized_Steve_v3_Kappa(paramsplus(1), paramsplus(2), kappa);
        %   growth_plus(jj, kk) = Kappa_Disp_Relation(paramsplus(1), paramsplus(2), init_guess_plus);
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
save(['Dispersion_Rate_Kappa_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'data_par.mat'])

%exit
