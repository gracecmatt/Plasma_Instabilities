%% 1D Lorentzian Solution Script Computing Relative Error in gamma
clear; clc; 
%%% --- parameters to change when tuning the numerical algorithm ---
% L (zetaph) = 1
% M (spectral) = 2^10
% Vmax (spectral) = 140

N = 40;
kARR = linspace(0.25,2,N); % k array
sigARR = linspace(0.5,4,N); % sigma array

errorXieIM = zeros(N,N); % relative error of imaginary part of root

mu = 10; % mean velocity

for i = 1:N % iterate through k array,
    k = kARR(i);

    for j = 1:N % iterate through sigma array
        sigma = sigARR(j);
    
        exact = mu*k+1 + 1i.*(-sigma*k);
    
        init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma,0) + mu*k;
        xi_guess_shiftscaled = (init_guess/k-mu)/sigma;
        omegaXieShiftScaled =  Lorentzian_Disp_Using_Xie(k*sigma,1,0,xi_guess_shiftscaled)*sigma*k + mu*k;
    
        errorXieIM(i,j) = abs(imag(exact)-imag(omegaXieShiftScaled))./abs(imag(exact));
        if errorXieIM(i,j) == 0; errorXieIM(i,j) = 10^(-16); end
    end
end

%%
mean(errorXieIM,"all")
std(errorXieIM,0,"all")
max(errorXieIM,[],"all")