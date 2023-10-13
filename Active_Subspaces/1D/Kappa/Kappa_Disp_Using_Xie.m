function gamma = Kappa_Disp_Using_Xie(k, sigma, mu, kappa, init_guess)

    % F = ['(gamma(',num2str(kappa),'+1)/gamma(',num2str(kappa),'+0.5)/sqrt(pi*(',...
    %     num2str(kappa),'-0.5)*',num2str(theta),'^2))*(1+(v-',num2str(mu),').^2/(',num2str(theta),...
    %     '^2*(',num2str(kappa),'-0.5))).^(-',num2str(kappa),'-1)'];

    C1 = (pi*sigma^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5))*(pi*sigma^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
    F = [num2str(C1,16),'*(1+(v-',num2str(mu,16),').^2/((',num2str(kappa,16),'-1.5)*',num2str(sigma,16),'^2)).^(-',num2str(kappa,16),')'];

    options = optimoptions('fsolve','Display','off','Algorithm','trust-region','FiniteDifferenceType','central');
    % options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
    Fn = 0; % option to define your own F
    N = 256; % number of Fourier coefficients to take
    D = @(omega) 1-1/k^2*zetaph(omega/k, Fn, F, N); % dielectric function
    omega = fsolve(D, init_guess, options);
    gamma = imag(omega);

end