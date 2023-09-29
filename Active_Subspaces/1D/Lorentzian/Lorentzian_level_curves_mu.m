clear; clc; close all;
mu0 = 0; muf = 900;
muplot = mu0:(muf-mu0)/499:muf;
count = 1;
initial_guesses = zeros(1,length(muplot));
omega_xie_rescaled = zeros(1,length(muplot));
omega_xie = zeros(1,length(muplot));

sigma = 1;
k = 0.5;

for mu=muplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0);
    initial_guesses(count) = init_guess + mu*k;
    omega_xie(count) = Maxwellian_Disp_Using_Xie(k, sigma, mu, init_guess+mu*k)*k;
    omega_xie_rescaled(count) = Maxwellian_Disp_Using_Xie(k*sigma, 1, 0, init_guess)*k*sigma + mu*k;
    count = count+1;
end

%% Figures
txt = ['$\sigma$ = ',num2str(sigma),', $k$ = ', num2str(k)];

figure
plot(muplot, imag(initial_guesses),'o-'); hold on
plot(muplot, imag(omega_xie),'.-');
plot(muplot, imag(omega_xie_rescaled),'linewidth',2);
title('Maxwellian - $\mu$ vs. $\gamma$','Interpreter','latex','FontSize',16)
xlabel('$\mu$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\mu)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(muplot(2)-muplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',14)

figure
plot(muplot, real(initial_guesses),'o-'); hold on
plot(muplot, real(omega_xie),'.-');
plot(muplot, real(omega_xie_rescaled),'linewidth',2);
title('Maxwellian - $\mu$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$\mu$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\mu)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','Location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(muplot(2)-muplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',14)