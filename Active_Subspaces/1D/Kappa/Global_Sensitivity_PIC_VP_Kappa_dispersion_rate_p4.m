close all;
rng('shuffle');
parpool(8);
% Initialize algorithm parameters
N = 512;                              %Number of samples for each parameter
h = 1e-6;                                      %Finite difference step size
kappa = 2;

% Pre-allocate memory
Nparams = 4;
growth = zeros(N,1);                      %Output of interest (growth rate)
growth_plus = zeros(N,Nparams);               %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
Xs = zeros(N,Nparams);                   %To save the normalized parameters
omega_error = zeros(N,1);           %Compare Xie/dielectric funtion results

% vals = [k, sigma, mu, kappa]
setvals = [0.5; 2; 1; kappa];

var = 0.05; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fixing mu
if setvals(3)==0
    xl(3) = -var;
    xu(3) = var;
end
if kappa==2 % keep kappa in the far-equilibrium region (Livadiotis & McComas 2013)
    xu(4) = min(xu(4),2.49); % keep kappa < 2.5
    xl(4) = max(xl(4),1.51); % keep kappa > 1.5
elseif kappa==6 % keep kappa in the near-equilibrium region
    xl(4) = max(xl(4),2.51); % keep kappa > 2.5
end

% Run simulation
tic
rng(sum(100*clock));
Xs = 2*rand(N^2,Nparams) - 1; % do sampling in serial
parfor jj = 1:N
    % Randomly sample parameters within acceptable ranges
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
    
    % Numerically solve 1D Vlasov-Poisson with randomly drawn parameters
    % init_guess = BohmGross_Kap(params(1),params(2),0,kappa);
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(params(1),params(2),0,params(4));
    xi_guess = init_guess/(params(1)*params(2)); % shifted and scaled

    omega = Kappa_Disp_Using_Xie(params(1)*params(2),1,0,params(4),xi_guess)*params(1)*params(2) + params(3)*params(1);
    dielectric = Kappa_dielectric(params(1),params(2),0,kappa,init_guess) + params(3)*params(1);
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
        
        % init_guess_plus = BohmGross_Kap(paramsplus(1),paramsplus(2),0,kappa);
        init_guess_plus = Vlasov_1D_linearized_Steve_v4_Kappa(paramsplus(1),paramsplus(2),0,paramsplus(4));
        xi_guess_plus = init_guess_plus/(paramsplus(1)*paramsplus(2)); % shifted and scaled
    
        omega0_plus = Kappa_Disp_Using_Xie(paramsplus(1)*paramsplus(2),1,0,paramsplus(4),xi_guess_plus)*paramsplus(1)*paramsplus(2);
        growth_plus(jj, kk) = imag(omega0_plus);
    end
end
parfor jj = 1:N
    % Calculate the appx gradients using finite differences
    grad_growth(:,jj) = (growth_plus(jj, :) - growth(jj))/h;
end
toc;

% Compute the singular value decomposition of the matrix C
[U,S,V] = svd(1/sqrt(N)*grad_growth);
w = U(:,1);
w2 = U(:,2);
    
%Compute the eigenvalues of C
evalues = diag(S.^2);
    
% Compute the condition number (need a new name??)
cond = evalues(1)/sum(evalues);

% Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));
   
%Save the trial data
save(['Dispersion_Kappa' int2str(kappa) '_P' int2str(Nparams) '_N' int2str(N) '_var' int2str(var*100) '_data.mat'])

%exit
delete(gcp('nocreate'))