kplot = 0.1:0.01:1;
count=1;
% gammas_root_finding = zeros(1,length(kplot));
initial_guesses = zeros(1,length(kplot));
gammas_xie = zeros(1,length(kplot));

kappa = 1;
theta = [1, 1];
mu = [0, 4];
beta = 0.9;

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k, theta, mu, beta, kappa);
    initial_guesses(count) = init_guess;
    % gammas_root_finding(count) = dielectric(k, theta(1), theta(2), mu(1), mu(2), beta, init_guess);
    gammas_xie(count) = Kappa_Bump_Disp_Using_Xie(k, theta(1), theta(2), mu(1), mu(2), beta, kappa, init_guess);
    count = count+1;
end

figure
plot(kplot, imag(initial_guesses))
hold on
% plot(kplot, gammas_root_finding)
plot(kplot, imag(gammas_xie))
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
title('Kappa with Bump - $k$ vs. $\gamma$ ($\kappa=1$)','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Root Finding with Xie Code')

figure
hold on
plot(kplot, real(initial_guesses));
plot(kplot, real(gammas_xie));
title('Kappa with Bump - $k$ vs. $\Omega$ ($\kappa=1$)','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')