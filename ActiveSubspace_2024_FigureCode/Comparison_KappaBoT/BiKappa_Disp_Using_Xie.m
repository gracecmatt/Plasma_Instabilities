function xi = BiKappa_Disp_Using_Xie(k, sigma1, sigma2, mu, v0, beta, kappa, init_guess)

    % =================== options for root finding ========================
    % options = optimoptions('fsolve','Display','none','Algorithm','trust-region','FiniteDifferenceType','central');
    options = optimoptions('fsolve','Display','off');

    % ========= compute gamma using Fourier series approximation ==========
    % kappa distribution has bump at mu-v0 and mu+v0
    F = [num2str(beta,16),'*(gamma(',num2str(kappa,16),')/gamma(',num2str(kappa,16),'-1/2)/sqrt(pi*(',...
            num2str(kappa,16),'-3/2)*',num2str(sigma1,16),'^2))*(1+(v-',num2str(mu,16),').^2/(',...
            num2str(sigma1,16),'^2*(',num2str(kappa,16),'-3/2))).^(-',num2str(kappa,16),') + ',...
        num2str(1-beta,16),'*(gamma(',num2str(kappa,16),')/gamma(',num2str(kappa,16),'-1/2)/sqrt(pi*(',...
            num2str(kappa,16),'-3/2)*',num2str(sigma2,16),'^2))*( (1+(v-',num2str(mu+v0,16),').^2/(',...
            num2str(sigma2,16),'^2*(',num2str(kappa,16),'-3/2))).^(-',num2str(kappa,16),') + (1+(v-',num2str(mu-v0,16),').^2/(',...
            num2str(sigma2,16),'^2*(',num2str(kappa,16),'-3/2))).^(-',num2str(kappa,16),') )/2'];
    Fn = 0; % 0 = option to define your own F
    Nfourier = 1000; % number of Fourier coefficients to take
    % D = @(z) 1-1/k^2*zetaph(z, Fn, F, Nfourier); % dielectric function
    D = @(z) k^2-zetaph(z, Fn, F, Nfourier); % dielectric function
    
    % ============================== solve ================================
    xi = fsolve(D, init_guess, options);
end