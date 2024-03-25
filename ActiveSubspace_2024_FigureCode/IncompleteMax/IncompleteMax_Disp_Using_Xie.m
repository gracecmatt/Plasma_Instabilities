function xi = IncompleteMax_Disp_Using_Xie(k, sigma, mu, nu, init_guess)

    % =================== options for root finding ========================
    options = optimoptions('fsolve','Display','off');

    % ========= compute omega using Fourier series approximation ==========
    F = ['2.*heaviside(real(v)-',num2str(nu,16),')./sqrt(pi*',num2str(sigma,16),'^2).*exp(-(v-',num2str(mu,16),').^2/(',num2str(sigma,16),'^2))./(1+erf(',num2str(mu,16),'-',num2str(nu,16),')./',num2str(sigma,16),'))'];

    Fn = 5; % 5 = option for incomplete maxwellian
    Nfourier = 2000; % number of Fourier coefficients to take
    D = @(z) 1-1/k^2*zetaph_IMax(z, Fn, F, Nfourier, mu, sigma, nu); % dielectric function

    xi = fsolve(D, init_guess, options); % root finder

end