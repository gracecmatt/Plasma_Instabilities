rng('shuffle');
parpool(4);%parpool(24);%was 32
% Initialize algorithm parameters
N = 30; %512   %Number of samples for each parameter
h = 1e-6;    %Finite difference step size
%trial = 1;   %Trial number (used when saving figures)

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

% Set upper and lower bounds for parameters - beta_2 = 1 - beta_1
% vals =  [beta_1; sigma2_1; mu_1; sigma2_2; mu_2; k]; 
%setvals = [0.9; 0.5; 0.2; 0.5; 3.8; 0.5];
% vals =  [k; mu_1; mu_2; sigma2_1; sigma2_2; beta_1]; 
setvals = [0.5; 0; 4; 0.25; 0.25; 0.8];

var = 0.25; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

xl(2) = -var;
xu(2) = var;

xu(6) = min(xu(6),1);

% xl = [0.4; -0.1; 3; 0.1; 0.1; 0.1];
% xu = [0.6;  0.1; 5; 0.5; 0.5; 0.9];  

% Run simulation
tic
parfor jj = 1:N
    rng(sum(100*clock)+pi*jj);
    % Randomly sample parameters within acceptable ranges
    Xs(jj,:) = 2*rand(1,Nparams) - 1;
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
    % Numerically solve 1D Vlasov-Poisson with baseline parameters
    growth(jj) = dispersion_growthrate_BiMax(params);
    %while growth(jj) < 0 || growth(jj) > 2 
    while (growth(jj) > 5  || growth(jj) < 1e-10)
        Xs(jj,:) = 2*rand(1,Nparams) - 1;
        params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
        growth(jj) = dispersion_growthrate_BiMax(params);
    end 
end

parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);                     
        % Numerically solve 1D Vlasov-Poisson with perturbed parameters
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        growth_plus(jj, kk)= dispersion_growthrate_BiMax(paramsplus);
        %while growth_plus(jj,kk) <0 || growth_plus(jj,kk) >2 % set to 2 for 10%, set to 1 for 25% runs
        while (growth_plus(jj,kk) > 5 || growth_plus(jj,kk) < 1e-10)
            xplus = randparams + h*I(:,kk);
            paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
            growth_plus(jj, kk)= dispersion_growthrate_BiMax(paramsplus);
        end 
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
save(['Dispersion_Rate_BiMax_P' int2str(Nparams) '_N' int2str(N) '_' num2str(var) 'data_par.mat'])

%exit