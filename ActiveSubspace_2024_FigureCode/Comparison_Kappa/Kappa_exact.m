function omega = Kappa_exact(k,kappa,init_guess)

    % ========= compute Z, Zp using Gauss hypergeometric function ==========
    % x = @(z) 1/2*(1-z/(1i*sqrt(kappa)));
    % Zgauss = @(z) 1-1/k^2*1i*(kappa+1/2)*(kappa-1/2)/(kappa^(3/2)*(kappa+1))*hypergeom( [1, 2*kappa+2], kappa+2, x(z) );
    Zpgauss = @(z) -(kappa+1/2)*(kappa-1/2)/(kappa^2*(kappa+2))*hypergeom( [2, 2*kappa+3], kappa+3, 1/2*(1-z/(1i*sqrt(kappa))) ) ;
    D = @(z) 1-1/k^2*Zpgauss(z);

    options = optimoptions('fsolve','Display','off');
    omega = fsolve(D, init_guess, options);
end