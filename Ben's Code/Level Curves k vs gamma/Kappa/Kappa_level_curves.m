kplot = 0.1:0.01:1;
count=1;
initial_guesses = zeros(1,length(kplot));
gamma_xie = zeros(1,length(kplot));

kappa = 1;
theta = 1;
mu = 0;

tic;

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v3_Kappa(k, kappa, mu, theta);
    initial_guesses(count) = init_guess;
    gamma_xie(count) = Kappa_Disp_Using_Xie(k, mu, theta, kappa, init_guess);
    count = count+1;
end

toc

% Create Plots
figure
plot(kplot, imag(initial_guesses))
hold on
plot(kplot, imag(gamma_xie))
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
title('Kappa Distribution ($\kappa=1$)','Interpreter','latex')
legend('Spectral Method','Dispersion Relation with Xie Code')

figure
plot(kplot, real(initial_guesses));
hold on
plot(kplot, real(gamma_xie));
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
title('Real Part of Initial Guess and Xie')
legend('Spectral Method','Dispersion Relation with Xie Code')
