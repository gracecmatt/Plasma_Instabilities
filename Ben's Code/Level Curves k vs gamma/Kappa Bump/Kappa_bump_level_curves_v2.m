kplot = 1:-0.01:0.1;
gamma_spectral_method = zeros(1,length(kplot));  % Still compute spectral method to compare
gamma_xie = zeros(1,length(kplot));

kappa = 1;
theta = [1, 1];
mu = [0, 4];
beta = 0.9;

% First for k=1 to get accurate initial guess
gamma_spectral_method(1) = Vlasov_1D_linearized_Steve_v3_Kappa(kplot(1), theta, mu, beta, kappa);
gamma_xie(1) = Kappa_Bump_Disp_Using_Xie(kplot(1), theta(1), theta(2), mu(1), mu(2), beta, imag(gamma_spectral_method(1)));

for count=2:length(kplot)
    k = kplot(count);
    gamma_spectral_method(count) = Vlasov_1D_linearized_Steve_v3_Kappa(kplot(count), theta, mu, beta, kappa);
    initial_guess = real(gamma_xie(count-1)) + 1i*imag(gamma_spectral_method(count));
    gamma_xie(count) = Kappa_Bump_Disp_Using_Xie(kplot(count), theta(1), theta(2), mu(1), mu(2), beta, initial_guess);
end

figure
plot(kplot, imag(gamma_spectral_method))
hold on
plot(kplot, imag(gamma_xie))
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
title('Kappa with Bump - $k$ vs. $\gamma$','Interpreter','latex','FontSize',14)
legend('Spectral Method', 'Xie Code (with different initial guess)')

figure
plot(kplot, real(gamma_xie))