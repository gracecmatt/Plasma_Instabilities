clear; clc; close all
kplot = 0.5:0.01:1.5;
count = 1;
initial_guesses = zeros(1,length(kplot));
omega_xie = zeros(1,length(kplot));

sigma = 1;
mu = 0;
nu = -1;
M = 1;
% vplot = -2.5:0.005:2.5;
% f0 = (atan(M*(vplot-nu))/pi+1/2).*exp(-(vplot-mu).^2/sigma^2)/sqrt(pi*sigma^2);
% plot(vplot,f0)

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0, nu, M) + mu*k;
    initial_guesses(count) = init_guess;
    omega_xie(count) = IncompleteMax_Disp_Using_Xie(k, sigma, 0, nu, M, init_guess) + mu*k;
    count = count+1;
end

figure
plot(kplot, imag(initial_guesses)); hold on
plot(kplot, imag(omega_xie));
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
title('Kappa with Bump - $k$ vs. $\gamma$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')

figure
plot(kplot, real(initial_guesses)); hold on
plot(kplot, real(omega_xie));
title('Kappa with Bump - $k$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')