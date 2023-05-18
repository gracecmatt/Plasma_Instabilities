function gamma = BiMaxwellian_Disp_Using_Xie(k,sigma1,sigma2,mu1,mu2,beta,init_guess)

    beta

    F = [num2str(beta),'/sqrt(2*pi*',num2str(sigma1),'^2)*exp(-abs(v-',num2str(mu1),')^2/(2*',num2str(sigma1),'^2)) + ',...
        num2str(1-beta),'/sqrt(2*pi*',num2str(sigma2),'^2)*exp(-abs(v-',num2str(mu2),')^2/(2*',num2str(sigma2),'^2))'];


    options = optimoptions('fsolve','Display','none','Algorithm','trust-region','FiniteDifferenceType','central');
    % options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
    Fn = 4; % option to define your own F
    NF = 32; % number of Fourier coefficients to take
    D = @(omega) 1-1/k^2*zetaph(omega/k, Fn, F, NF); % dielectric function
    omega = fsolve(D, init_guess, options);
    gamma = imag(omega);

end