kplot = 0.1:0.005:1;
count=1;
gamma_spectral_method = zeros(1,length(kplot));
gamma_root_finding = zeros(1,length(kplot));
gamma_xie = zeros(1,length(kplot));

kappa = 1;
theta = 1;
mu = 0;

tic;

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v3_Kappa(k, kappa, theta);
    gamma_spectral_method(count) = imag(init_guess);
    gamma_root_finding(count) = dielectric(k, theta, theta, mu, mu, 1, init_guess);
    gamma_xie(count) = Kappa_Disp_Using_Xie(k, theta, mu, init_guess);
    count = count+1;
end

toc

% Create Plots
figure
plot(kplot, gamma_spectral_method)
hold on
plot(kplot, gamma_root_finding)
plot(kplot, gamma_xie)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
title('Kappa Distribution ($\kappa=1$)','Interpreter','latex')
legend('Spectral Method','Old Dispersion Relation','Dispersion Relation with Xie Code')

% Plotting Error
xie_error = max(gamma_xie - gamma_root_finding)
