clear; clc; close all
M0 = 1; Mf=70;
Mplot = M0:(Mf-M0)/99:Mf;
count = 1;
initial_guesses = zeros(1,length(Mplot));
omega_xie = zeros(1,length(Mplot));

sigma = 1;
mu = 0;
nu = -1;
k = 0.5;
% vplot = -2.5:0.005:2.5;
% f0 = (atan(M*(vplot-nu))/pi+1/2).*exp(-(vplot-mu).^2/sigma^2)/sqrt(pi*sigma^2);
% plot(vplot,f0)

for M=Mplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0, nu, M) + mu*k;
    initial_guesses(count) = init_guess;
    omega_xie(count) = IncompleteMax_Disp_Using_Xie(k, sigma, 0, nu, M, init_guess) + mu*k;
    count = count+1;
end
%% Figures

% Imaginary part
figure
subplot(1,2,1)
plot(Mplot, imag(initial_guesses)); hold on
plot(Mplot, imag(omega_xie));
xlabel('$M$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(M)$','Interpreter','latex','FontSize',16)
title('$M$ vs. $\gamma$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','location','Best')

% Real plot
% figure
subplot(1,2,2)
plot(Mplot, real(initial_guesses)); hold on
plot(Mplot, real(omega_xie));
title('$M$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$M$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(M)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','location','Best')