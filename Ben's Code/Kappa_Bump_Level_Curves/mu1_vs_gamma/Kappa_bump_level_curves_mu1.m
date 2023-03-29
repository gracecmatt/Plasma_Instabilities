mu1plot = -2:0.04:2;
count=1;
initial_guesses = zeros(1,length(mu1plot));
gammas_xie = zeros(1,length(mu1plot));
   
k = 0.5;
kappa = 1;
theta1 = 1;
theta2 = 1;
mu2 = 4;
beta = 0.9;

for mu1=mu1plot
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k, theta1, theta2, mu1, mu2, beta, kappa);
    initial_guesses(count) = init_guess;
    gammas_xie(count) = Kappa_Bump_Disp_Using_Xie(k, theta1, theta2, mu1, mu2, beta, kappa, init_guess);
    count = count+1;
end

figure
plot(mu1plot, imag(initial_guesses))
hold on
% plot(kplot, gammas_root_finding)
plot(mu1plot, imag(gammas_xie))
xlabel('$\mu_1$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\mu_1)$','Interpreter','latex','FontSize',16)
title('Kappa with Bump - $\mu_1$ vs. $\gamma$ ($\kappa=1$)','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Root Finding with Xie Code')

figure
hold on
plot(mu1plot, real(initial_guesses));
plot(mu1plot, real(gammas_xie));
title('Kappa with Bump - $\mu_1$ vs. $\Omega$ ($\kappa=1$)','Interpreter','latex','FontSize',16)
xlabel('$\mu_1$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\mu_1)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')