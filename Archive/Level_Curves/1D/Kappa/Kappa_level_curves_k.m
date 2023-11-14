kplot = 0.1:0.01:1;
count=1; % increment variable
initial_guesses = zeros(1,length(kplot));
gamma_xie = zeros(1,length(kplot));
gamma_root = zeros(1,length(kplot));

kappa = 1; % choose from {1,2,6}
theta = 1;
mu = 0;

tic;

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v3_Kappa(k, theta, mu, kappa); % spectral method solution
    initial_guesses(count) = init_guess;
    gamma_xie(count) = Kappa_Disp_Using_Xie(k, theta, mu, kappa, init_guess); % uses Fourier series approximation of integral
    gamma_root(count) = dielectric_kappa(k,theta,theta,mu,mu,1,init_guess); % integral computed in Mathematica
    count = count+1;
end

toc

% Create Plots
figure
plot(kplot, imag(initial_guesses))
hold on
plot(kplot, imag(gamma_xie))
plot(kplot, imag(gamma_root))
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
title(['Kappa Distribution ($\kappa=',num2str(kappa),'$)'],'Interpreter','latex')
legend('Spectral Method','Dispersion Relation with Xie Code','Root Finding with dielectric.m')

figure
plot(kplot, real(initial_guesses));
hold on
plot(kplot, real(gamma_xie));
plot(kplot, real(gamma_root));
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
title('Real Part of Initial Guess and Xie')
legend('Spectral Method','Dispersion Relation with Xie Code','Root Finding with dielectric.m')