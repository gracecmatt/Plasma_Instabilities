clear; clc; 
sigma0 = 1; sigmaf = 2;
sigmaplot = sigma0:(sigmaf-sigma0)/199:sigmaf;
count = 1;
initial_guesses = zeros(1,length(sigmaplot));
omega_xie_rescaled = zeros(1,length(sigmaplot));
omega_xie = zeros(1,length(sigmaplot));

mu = 100;
k = 0.5;
% sigma = sigmaf/sigmaf;
% vplot = [-fliplr(logspace(-5,0,200)),logspace(-5,0,199)]/sigmaf;
% f0 = exp(-(vplot-mu).^2/sigma^2)/sqrt(pi*sigma^2);
% semilogy(vplot,f0)

for sigma=sigmaplot
    % init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0);
    init_guess = analyticAPPROX(k,sigma,mu);
    initial_guesses(count) = init_guess;
    omega_xie(count) = Maxwellian_Disp_Using_Xie(k, sigma, mu, init_guess)*k;
    omega_xie_rescaled(count) = Maxwellian_Disp_Using_Xie(k*sigma, 1, 0, init_guess)*k*sigma + mu*k;
    count = count+1;
end

%% Figures
txt = ['$\mu$ = ',num2str(mu),', $k$ = ', num2str(k)];

figure
plot(sigmaplot, imag(initial_guesses),'o-'); hold on
plot(sigmaplot, imag(omega_xie),'.-'); hold on
plot(sigmaplot, imag(omega_xie_rescaled));
title('Maxwellian - $\sigma$ vs. $\gamma$','Interpreter','latex','FontSize',16)
xlabel('$\sigma$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\sigma)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(sigmaplot(2)-sigmaplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',14)

figure
plot(sigmaplot, real(initial_guesses),'o-'); hold on
plot(sigmaplot, real(omega_xie),'.-'); hold on
plot(sigmaplot, real(omega_xie_rescaled));
title('Maxwellian - $\sigma$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$\sigma$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\sigma)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','Location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(sigmaplot(2)-sigmaplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',14)