function xi = Kappa_Disp_Using_Xie(k, sigma, mu, kappa, init_guess)

    % =================== options for root finding ========================
    % options = optimoptions('fsolve','Display','none','Algorithm','trust-region','FiniteDifferenceType','central');
    options = optimoptions('fsolve','Display','off');

    % ========= compute gamma using Fourier series approximation ==========
    F = ['(gamma(',num2str(kappa,16),')/gamma(',num2str(kappa,16),'-0.5)/sqrt(pi*(',...
        num2str(kappa,16),'-1.5)*',num2str(sigma,16),'^2))*(1+(v-',num2str(mu,16),').^2/(',num2str(sigma,16),...
        '^2*(',num2str(kappa,16),'-1.5))).^(-',num2str(kappa,16),')'];

    Fn = 0; % 0 = option to define your own F
    Nfourier = 900; % number of Fourier coefficients to take
    D = @(z) 1-1/k^2*zetaph(z, Fn, F, Nfourier); % dielectric function
    
    xi = fsolve(D, init_guess, options);

end