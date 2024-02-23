rng('shuffle');
parpool(12);%parpool(24);%was 32

% Initialize algorithm parameters
N = 600;   %Number of samples for each parameter
h = 1e-6;  %Finite difference step size

% Pre-allocate memory
growth = zeros(N,1);                        %Output of interest (growth rate)
init_guesses = zeros(N,1);                %Output of interest (growth rate)
Nparams = 6;
growth_plus = zeros(N,Nparams);               %Perturbed output of interest
grad_growth = zeros(Nparams,N);             %Gradient of output of interest 
% Xs = zeros(N,Nparams);                    %To save the normalized paramters
                    
% Set upper and lower bounds for parameters
%vals =  [k; mu1; mu2; sigma1^2; sigma2^2; beta]; 
setvals = [0.5; 0; 4; 0.25; 0.25; 0.8];

var = 0.25; % x 100% variation considered 
xl = (1-var)*setvals;
xu = (1+var)*setvals;

% fix mu around 0
if setvals(2) == 0 
    xl(2) = -var;
    xu(2) = var;
end

xu(6) = min(xu(6),1);

% Run simulation
tic
rng(sum(100*clock)+pi); 
Xs = 2*rand(N,Nparams) - 1;
%Baseline Parameters
parfor jj = 1:N
    % Randomly sample parameters within acceptable ranges
    params = 1/2*(diag(xu - xl)*Xs(jj,:)' + (xu + xl));
    % Compute root of dispersion relation
    init_guess = Vlasov_1D_linearized_Steve_BiMax_v4(params);

    xi_guess = init_guess/(params(1)*sqrt(params(4))); % shifted and scaled
 
    omega = BiMaxwellian_Disp_Using_Xie(params(1)*sqrt(params(4)),1,params(5)/params(4),0,(params(3)-params(2))/sqrt(params(4)),params(6),xi_guess)*params(1)*sqrt(params(4));
    growth(jj) = imag(omega);
end

%Perturbed Parameters 
parfor jj = 1:N
    randparams = Xs(jj,:)';
    for kk = 1:Nparams
        I = eye(Nparams);                     
        % Compute root of dispersion relation
        xplus = randparams + h*I(:,kk);
        paramsplus = 1/2*(diag(xu - xl)*xplus + (xu + xl));
        init_guess_plus = Vlasov_1D_linearized_Steve_BiMax_v4(paramsplus);

        % xi_guess_plus = init_guess_plus/paramsplus(1); % shifted (not scaled)
        xi_guess_plus = init_guess_plus/(paramsplus(1)*sqrt(paramsplus(4))); % shifted and scaled
     
        % omega_plus = BiMaxwellian_Disp_Using_Xie(paramsplus(1),paramsplus(4),paramsplus(5),0,paramsplus(3)-paramsplus(2),paramsplus(6),xi_guess_plus)*paramsplus(1);
        omega_plus = BiMaxwellian_Disp_Using_Xie(paramsplus(1)*sqrt(paramsplus(4)),1,paramsplus(5)/paramsplus(4),0,(paramsplus(3)-paramsplus(2))/sqrt(paramsplus(4)),paramsplus(6),xi_guess_plus)*paramsplus(1)*sqrt(paramsplus(4));
        growth_plus(jj, kk) = imag(omega_plus);
    end
end
parfor jj = 1:N
    % Calculate the appx gradients using finite differences
    grad_growth(:,jj) = (growth_plus(jj, :) - growth(jj))/h;
        
end
toc

% Compute the weights and eigenvalues
%Compute the singular value decomposition of the matrix C
[U,S,V] = svd(1/sqrt(N)*grad_growth);
w = U(:,1);
w2 = U(:,2);
    
%Compute the eigenvalues of C
evalues = diag(S.^2);

% Compute the first two eigenvalue ratios
eta(1) = evalues(1)/sum(evalues);
eta(2) = (evalues(1)+evalues(2))/sum(evalues);
    
%Find the difference of max and min grad_growth to check for errors
diff_growth = max(max(grad_growth)) - min(min(grad_growth));
   
%Save the trial data
save(['Data/Dispersion_BiMax_P' int2str(Nparams) '_N' int2str(N) '_var' int2str(var*100) '_data.mat'])

%exit
delete(gcp('nocreate'))