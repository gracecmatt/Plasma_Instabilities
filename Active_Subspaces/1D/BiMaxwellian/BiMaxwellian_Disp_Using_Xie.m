function xi = BiMaxwellian_Disp_Using_Xie(k, sigma1, sigma2, mu1, mu2, beta, init_guess)

    F = [num2str(beta,16),'/sqrt(pi*',num2str(sigma1,16),'^2)*exp(-(v-',num2str(mu1,16),').^2/(',num2str(sigma1,16),'^2)) + ',...
        num2str(1-beta,16),'/sqrt(pi*',num2str(sigma2,16),'^2)*exp(-(v-',num2str(mu2,16),').^2/(',num2str(sigma2,16),'^2))'];

    options = optimoptions('fsolve','Display','none');%,'Algorithm','trust-region','FiniteDifferenceType','central');
    % options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
    
    Fn = 0; % 0 = option to define your own F
    Nfourier = 1200; % number of Fourier coefficients to take
    D = @(z) 1-1/k^2*zetaph(z, Fn, F, Nfourier); % dielectric function
    
    xi = fsolve(D, init_guess, options); % root finder
    xiReal = real(xi);
    xiImag = imag(xi);

    xi = abs(xiReal)+1i*xiImag;

end