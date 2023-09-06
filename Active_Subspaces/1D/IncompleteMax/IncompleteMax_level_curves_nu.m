clear; clc; close all
nuplot = -2:0.01:2;
count = 1;
initial_guesses = zeros(1,length(nuplot));
omega_xie = zeros(1,length(nuplot));

sigma = 1;
mu = 0;
M = 1;
k = 1;
% vplot = -2.5:0.005:2.5;
% f0 = (atan(M*(vplot-nu))/pi+1/2).*exp(-(vplot-mu).^2/sigma^2)/sqrt(pi*sigma^2);
% plot(vplot,f0)

for nu=nuplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0, nu, M) + mu*k;
    initial_guesses(count) = init_guess;
    omega_xie(count) = IncompleteMax_Disp_Using_Xie(k, sigma, 0, nu, M, init_guess) + mu*k;
    count = count+1;
end

% Imaginary part
figure
plot(nuplot, imag(initial_guesses)); hold on
plot(nuplot, imag(omega_xie));
xlabel('$\nu$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\nu)$','Interpreter','latex','FontSize',16)
title('Incomplete Maxwellian - $\nu$ vs. $\gamma$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')

% Real plot
figure
plot(nuplot, real(initial_guesses)); hold on
plot(nuplot, real(omega_xie));
title('Incomplete Maxwellian - $\nu$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$\nu$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\nu)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding')