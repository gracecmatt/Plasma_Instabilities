function xi = BiMaxwellian_Disp_Using_Xie(k, sigma1, sigma2, mu1, mu2, beta, init_guess)
    % =================== options for root finding ========================
    % options = optimoptions('fsolve','Display','none','Algorithm','trust-region','FiniteDifferenceType','central');
    options = optimoptions('fsolve','Display','off');

    % ========= compute gamma using Fourier series approximation ==========
    F = [num2str(beta,16),'/sqrt(pi*',num2str(sigma1,16),'^2)*exp(-(v-',num2str(mu1,16),').^2/(',num2str(sigma1,16),'^2)) + ',...
        num2str(1-beta,16),'/sqrt(pi*',num2str(sigma2,16),'^2)*exp(-(v-',num2str(mu2,16),').^2/(',num2str(sigma2,16),'^2))'];
    Fn = 0; % 0 = option to define your own F
    Nfourier = 500; % number of Fourier coefficients to take
    D = @(z) 1-1/k^2*zetaph(z, Fn, F, Nfourier); % dielectric function

    % ============================== solve ================================
    %% Just use initial guess as is
    xi = fsolve(D, init_guess, options);

    %% Use +- real part of initial condition and pick root with largest gamma
    % init_guess_1 = abs(real(init_guess))+1i*imag(init_guess);
    % init_guess_2 = -abs(real(init_guess))+1i*imag(init_guess);
    % 
    % xi = fsolve(D, init_guess_1, options); % root finder with positive real part of initial guess
    % xiReal_1 = real(xi);
    % xiImag_1 = imag(xi);
    % 
    % xi = fsolve(D, init_guess_2, options); % root finder with negative real part of initial guess
    % xiReal_2 = real(xi);
    % xiImag_2 = imag(xi);
    % 
    % % choose the solution with the largest "gamma" or imaginary part of the root
    % if(xiImag_1>xiImag_2) 
    %     xi = abs(xiReal_1)+1i*xiImag_1; 
    % elseif(xiImag_1<=xiImag_2)
    %     xi = -abs(xiReal_2)+1i*xiImag_2; 
    % end

end