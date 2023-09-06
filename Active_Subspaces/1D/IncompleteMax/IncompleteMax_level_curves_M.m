clear; clc; close all
Mplot = 1:50:10001;
count = 1;
initial_guesses = zeros(1,length(Mplot));
omega_xie = zeros(1,length(Mplot));

sigma = 1;
mu = 0;
nu = -1;
k = 1;
% vplot = -2.5:0.005:2.5;
% f0 = (atan(M*(vplot-nu))/pi+1/2).*exp(-(vplot-mu).^2/sigma^2)/sqrt(pi*sigma^2);
% plot(vplot,f0)

for M=Mplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0, nu, M) + mu*k;
    initial_guesses(count) = init_guess;
    omega_xie(count) = IncompleteMax_Disp_Using_Xie(k, sigma, 0, nu, M, init_guess) + mu*k;
    count = count+1;
end

% Imaginary part
figure
plot(Mplot, imag(initial_guesses)); hold on
plot(Mplot, imag(omega_xie));
xlabel('$M$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(M)$','Interpreter','latex','FontSize',16)
title('Incomplete Maxwellian - $M$ vs. $\gamma$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')

% Real plot
figure
plot(Mplot, real(initial_guesses)); hold on
plot(Mplot, real(omega_xie));
title('Incomplete Maxwellian - $M$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$M$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(M)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')