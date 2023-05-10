muplot = -2:0.02:2;
count=1;
initial_guesses = zeros(1,length(kplot));
gamma_xie = zeros(1,length(kplot));
gamma_root = zeros(1,length(kplot));

kappa = 1; % choose from {1,2,6}
theta = 1;
k = 0.5;

tic;

for mu=muplot
    init_guess = Vlasov_1D_linearized_Steve_v3_Kappa(k, theta, mu, kappa);
    initial_guesses(count) = init_guess;
    gamma_xie(count) = Kappa_Disp_Using_Xie(k, theta, mu, kappa, init_guess);
    gamma_root(count) = dielectric_kappa(k,theta,theta,mu,mu,1,init_guess);
    count = count+1;
end

toc

% Create Plots
figure
plot(muplot, imag(initial_guesses))
hold on
plot(muplot, imag(gamma_xie))
plot(muplot, imag(gamma_root))
xlabel('$\mu$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\mu)$','Interpreter','latex','FontSize',16)
title(['Kappa Distribution ($\kappa=',num2str(kappa),'$)'],'Interpreter','latex')
legend('Spectral Method','Dispersion Relation with Xie Code','Root Finding with dielectric.m')

figure
plot(muplot, real(initial_guesses));
hold on
plot(muplot, real(gamma_xie));
plot(muplot, real(gamma_root));
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\mu)$','Interpreter','latex','FontSize',16)
title('Real Part of Initial Guess and Xie')
legend('Spectral Method','Dispersion Relation with Xie Code','Root Finding with dielectric.m')
