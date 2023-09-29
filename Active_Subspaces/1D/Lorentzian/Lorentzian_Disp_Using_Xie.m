function xi = Lorentzian_Disp_Using_Xie(k, sigma, mu, init_guess)

    % =================== options for root finding ========================
    options = optimoptions('fsolve','Display','none');%,'Algorithm','trust-region','FiniteDifferenceType','central');
    % options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');


    % ========= compute gamma using Fourier series approximation ==========
    % F = ['1/sqrt(pi*',num2str(sigma,16),'^2)*exp(-(v-',num2str(mu,16),').^2/(',num2str(sigma,16),'^2))'];
    F = ['1/pi*',num2str(sigma,16),'./((v-',num2str(mu,16),').^2+',num2str(sigma,16),'^2)'];
    Fn = 0; % 0 = option to define your own F
    Nfourier = 1500; % number of Fourier coefficients to take
    D = @(z) 1-1/k^2*zetaph(z, Fn, F, Nfourier); % dielectric function

    xi = fsolve(D, init_guess, options); % root finder
    xiReal = real(xi);
    xiImag = imag(xi);
    xi = abs(xiReal) + 1i*xiImag;


    % % ============== compute gamma using plasma Z function ================
    % % % sigma = standard deviation
    % A = @(u) (1/(sqrt(2)*sigma))*(u - mu);
    % Z = @(u) -(1/sigma^2)*(1+A(u)*zetaf(A(u)));
    % % Z_np = @(u) 1/sqrt(2*sigma^2)*zetaf(A(u)); % formula with "no prime" in equilibrium distribution
    % D_zeta = @(omega) 1-1/k^2*Z(omega); % dielectric function
    % 
    % omega_zetaf = fsolve(D_zeta, init_guess, options);
    % gamma_zetaf = k*imag(omega_zetaf);


    % ================ determine error between methods ====================
    % diff_Zp = abs(Z_np(omega_zetaf) - zetaph(omega_zetaf, Fn, F, Nfourier))/abs(Z_np(omega_zetaf));
    % diff_gamma = abs((gamma_zetaf - gamma_fourier)/gamma_zetaf);
    % g=2;

    % gamma = gamma_fourier;

end