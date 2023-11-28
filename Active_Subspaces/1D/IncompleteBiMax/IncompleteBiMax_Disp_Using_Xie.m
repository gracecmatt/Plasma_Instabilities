function xi = IncompleteBiMax_Disp_Using_Xie(k, sigma1, sigma2, mu1, mu2, beta, nu, init_guess)

    % =================== options for root finding ========================
    options = optimoptions('fsolve','Display','off');

    % ========= compute omega using Fourier series approximation ==========
    % 2*heaviside(v-nu)*( beta*exp(-(v-mu1)^2/sigma1^2)/sqrt(pi*sigma1^2)/(1+erf((mu1-nu)/sigma1)) 
    %               + (1-beta)*exp(-(v-mu2)^2/sigma2^2)/sqrt(pi*sigma2^2)/(1+erf((mu2-nu)/sigma2)) )
    F = ['2.*heaviside(real(v)-',num2str(nu,16),').*( ',num2str(beta,16),'.*exp(-(v-',num2str(mu1,16),').^2./',num2str(sigma1,16),'^2)./sqrt(pi*',num2str(sigma1,16),'^2)./(1+erf((',num2str(mu1,16),'-',num2str(nu,16),')/',num2str(sigma1,16),'))'... 
                                                 '+ ',num2str(1-beta,16),'.*exp(-(v-',num2str(mu2,16),').^2./',num2str(sigma2,16),'^2)./sqrt(pi*',num2str(sigma2,16),'^2)./(1+erf((',num2str(mu2,16),'-',num2str(nu,16),')/',num2str(sigma2,16),')) )'];

    Fn = 5; % 5 = option for incomplete maxwellian
    Nfourier = 2000; % number of Fourier coefficients to take
    D = @(z) 1-1/k^2*zetaph_IBiMax(z, Fn, F, Nfourier, sigma1, sigma2, mu1, mu2, beta, nu); % dielectric function

    xi = fsolve(D, init_guess, options); % root finder

end