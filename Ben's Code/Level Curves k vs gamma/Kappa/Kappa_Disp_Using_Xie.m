function gamma = Kappa_Disp_Using_Xie(k, theta, mu, init_guess)

kappa = 1; % fixed kappa
F = ['(gamma(',num2str(kappa),'+1)/gamma(',num2str(kappa),'+0.5)/sqrt(pi*(',...
        num2str(kappa),'-0.5)*',num2str(theta),'^2))*(1+(v-',num2str(mu),').^2/(',num2str(theta),...
        '^2*(',num2str(kappa),'-0.5))).^(-',num2str(kappa),'-1)'];

options = optimoptions('fsolve','Display','off');
z = fsolve(@(omega) 1-1/k^2*zetaph(omega/k, 0, F, 256), init_guess, options);
gamma = imag(z);

end