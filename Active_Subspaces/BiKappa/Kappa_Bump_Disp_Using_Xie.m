function gamma = Kappa_Bump_Disp_Using_Xie(k, theta1, theta2, mu1, mu2, beta, kappa, init_guess)

F = [num2str(beta),'*((gamma(',num2str(kappa),'+1)/gamma(',num2str(kappa),'+0.5)/sqrt(pi*(',...
        num2str(kappa),'-0.5)*',num2str(theta1),'^2))*(1+(v-',num2str(mu1),').^2/(',num2str(theta1),'^2*(',...
        num2str(kappa),'-0.5))).^(-',num2str(kappa),'-1)) + ',...
     num2str(1-beta),'*((gamma(',num2str(kappa),'+1)/gamma(',num2str(kappa),'+0.5)/sqrt(pi*(',...
        num2str(kappa),'-0.5)*',num2str(theta2),'^2))*(1+(v-',num2str(mu2),').^2/(',num2str(theta2),'^2*(',...
        num2str(kappa),'-0.5))).^(-',num2str(kappa),'-1))'];

% disp(k)
    options = optimoptions('fsolve','Display','none','Algorithm','trust-region','FiniteDifferenceType','central');
    % options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
    Fn = 0; % option to define your own F
    N = 512; % number of Fourier coefficients to take
    % z = omega/k; %root to solve for 
    % z = -0.5:0.01:0.5; % I really don't know what z is *****************************
    % [x,y]=meshgrid(-1:0.1:1,-1:0.1:1); z=x+1i*y; %[Zp,Z]=zetaph(z,1)
    % [Zp,~] = zetaph(z-omega/k, Fn, F, N);
    D = @(omega) 1-1/k^2*zetaph(omega/k, Fn, F, N);
    omega = fsolve(D, init_guess, options);
    gamma = imag(omega);

%    omega = fsolve(F, init_guess, options);


end