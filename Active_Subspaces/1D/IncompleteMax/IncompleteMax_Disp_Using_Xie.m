function omega = IncompleteMax_Disp_Using_Xie(k, sigma, mu, nu, M, init_guess)

    % =================== options for root finding ========================
    options = optimoptions('fsolve','Display','off');%,'Algorithm','trust-region','FiniteDifferenceType','central');
    % options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');

    % ========= compute gamma using Fourier series approximation ==========
    % F = ['heaviside(v-',num2str(nu,16),')./sqrt(pi*',num2str(sigma,16),'^2).*exp(-(v-',num2str(mu,16),').^2/(',num2str(sigma,16),'^2))'];
    F = ['(atan(',num2str(M,16),'.*(v-',num2str(nu,16),'))./pi+1/2).*exp(-(v-',num2str(mu,16),').^2/(',num2str(sigma,16),'^2))./sqrt(pi*',num2str(sigma,16),'^2)'];
    
    Fn = 0; % 0 = option to define your own F
    Nfourier = 1000; % number of Fourier coefficients to take
    D = @(omega) 1-1/k^2*zetaph(omega, Fn, F, Nfourier); % dielectric function

    omega = k*fsolve(D, init_guess, options); % root finder
    gamma = imag(omega);

end