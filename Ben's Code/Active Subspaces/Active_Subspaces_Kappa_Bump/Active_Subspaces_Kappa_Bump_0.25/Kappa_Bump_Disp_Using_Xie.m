function gamma = Kappa_Bump_Disp_Using_Xie(k, theta1, theta2, mu1, mu2, beta, kappa, init_guess)

F = [num2str(beta),'*((gamma(',num2str(kappa),'+1)/gamma(',num2str(kappa),'+0.5)/sqrt(pi*(',...
        num2str(kappa),'-0.5)*',num2str(theta1),'^2))*(1+(v-',num2str(mu1),').^2/(',num2str(theta1),'^2*(',...
        num2str(kappa),'-0.5))).^(-',num2str(kappa),'-1)) + ',...
     num2str(1-beta),'*((gamma(',num2str(kappa),'+1)/gamma(',num2str(kappa),'+0.5)/sqrt(pi*(',...
        num2str(kappa),'-0.5)*',num2str(theta2),'^2))*(1+(v-',num2str(mu2),').^2/(',num2str(theta2),'^2*(',...
        num2str(kappa),'-0.5))).^(-',num2str(kappa),'-1))'];

% disp(k)
options = optimoptions('fsolve','Display','none','Algorithm','trust-region','FiniteDifferenceType','central');
z = fsolve(@(omega) 1-1/k^2*zetaph(omega/k, 0, F, 512), init_guess, options);
gamma = imag(z);

end