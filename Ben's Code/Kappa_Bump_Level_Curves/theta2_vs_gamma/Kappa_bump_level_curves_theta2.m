theta2plot = 0.5:0.01:1.5;
count=1;
initial_guesses = zeros(1,length(theta1plot));
gammas_xie = zeros(1,length(theta1plot));

k = 0.5;
kappa = 1;
theta1 = 1;
mu1 = 0;
mu2 = 4;
beta = 0.9;

for theta2=theta2plot
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k, theta1, theta2, mu1, mu2, beta, kappa);
    initial_guesses(count) = init_guess;
    gammas_xie(count) = Kappa_Bump_Disp_Using_Xie(k, theta1, theta2, mu1, mu2, beta, kappa, init_guess);
    count = count+1;
end

figure
plot(theta2plot, imag(initial_guesses))
hold on
plot(theta2plot, imag(gammas_xie))
xlabel('$\theta_2$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\theta_2)$','Interpreter','latex','FontSize',16)
title('Kappa with Bump - $\theta_2$ vs. $\gamma$ ($\kappa=1$)','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Root Finding with Xie Code')

figure
hold on
plot(mu1plot, real(initial_guesses));
plot(mu1plot, real(gammas_xie));
title('Kappa with Bump - $\theta_2$ vs. $\Omega$ ($\kappa=1$)','Interpreter','latex','FontSize',16)
xlabel('$\theta_2$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\theta_2)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')