clear; clc
theta1plot = 0.001:0.001:0.01;
count = 1;
initial_guesses = zeros(1,length(theta1plot));
gammas_xie = zeros(1,length(theta1plot));
% gammas_exact = zeros(1,length(theta1plot));

k = 0.9;
kappa = 1.9; % kappa > 3/2
theta2 = 1;
mu1 = 800;
mu2 = 1;
beta = 1;
rmcount = [];

for theta1=theta1plot
    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k, theta1, theta2, mu1, mu2, beta, kappa);
    initial_guesses(count) = init_guess;
    % gammas_exact(count) = dielectricBOT_newkappa2(k, theta1, theta2, mu1, mu2, beta, kappa, init_guess);
    gammas_xie(count) = Kappa_Bump_Disp_Using_Xie(k, theta1, theta2, mu1, mu2, beta, kappa, init_guess);
    % if (abs(gammas_xie(count)) > 10) % || (abs(gammas_exact(count)) > 10)
    %    rmcount = [rmcount,count];
    % end
    count = count+1;
end

initial_guesses(rmcount) = [];
% gammas_exact(rmcount) = [];
gammas_xie(rmcount) = [];
theta1plot(rmcount) = [];

figure
plot(theta1plot, imag(initial_guesses),'.-','linewidth',1); hold on
plot(theta1plot, imag(gammas_xie),'.-','linewidth',1)
% plot(theta1plot, imag(gammas_exact),'.-','linewidth',1)
xlabel('$\theta_1$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\theta_1)$','Interpreter','latex','FontSize',16)
title('Kappa with Bump - $\theta_1$ vs. $\gamma$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding');%,'Exact Dielectric Function','location','best')

figure
plot(theta1plot, real(initial_guesses),'.-','linewidth',1); hold on
plot(theta1plot, real(gammas_xie),'.-','linewidth',1);
% plot(theta1plot, real(gammas_exact),'.-','linewidth',1);
title('Kappa with Bump - $\theta_1$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$\theta_1$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\theta_1)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding');%,'Exact Dielectric Function','location','best')