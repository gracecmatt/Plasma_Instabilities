function gamma = Kappa_Disp_Using_Xie(k, theta, mu, kappa, init_guess)

    F = ['(gamma(',num2str(kappa),'+1)/gamma(',num2str(kappa),'+0.5)/sqrt(pi*(',...
        num2str(kappa),'-0.5)*',num2str(theta),'^2))*(1+(v-',num2str(mu),').^2/(',num2str(theta),...
        '^2*(',num2str(kappa),'-0.5))).^(-',num2str(kappa),'-1)'];

    options = optimoptions('fsolve','Display','off','Algorithm','trust-region','FiniteDifferenceType','central');
    % options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
    Fn = 0; % option to define your own F
    N = 256; % number of Fourier coefficients to take
    D = @(omega) 1-1/k^2*zetaph(omega/k, Fn, F, N); % dielectric function
    omega = fsolve(D, init_guess, options);
    gamma = imag(omega);

end