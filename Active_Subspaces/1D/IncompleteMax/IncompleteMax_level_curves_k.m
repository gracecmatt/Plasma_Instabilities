clear; clc; close all
k0 = 0.57; kf = 1.2;
kplot = k0:(kf-k0)/99:kf;
count = 1;
omega_spectral = zeros(1,length(kplot));
omega_xie = zeros(1,length(kplot));
omega_xie_rescaled=zeros(1,length(kplot));

sigma = 1;
mu = 0;
nu = -0.9;
M = 10;
v = -2.5:0.005:2.5;
f0 = @(v) (1+erf(M*(v-nu))).*exp(-(v-mu).^2/sigma^2)/(sqrt(pi*sigma^2)*(1+erf(M*(mu-nu)/sqrt(1+M^2*sigma^2))));
% f0 = @(v) (atan(M*(v-nu))/pi+1/2).*exp(-(v-mu).^2/sigma^2)/sqrt(pi*sigma^2);
% figure; plot(v,f0(v))

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0, nu, M);
    xi_guess = (init_guess+mu*k)/k;
    xi_guess_rescaled = init_guess/(sigma*k);

    omega_spectral(count) = init_guess + mu*k;
    omega_xie(count) = IncompleteMax_Disp_Using_Xie(k, sigma, mu, nu, xi_guess)*k;
    omega_xie_rescaled(count)= IncompleteMax_Disp_Using_Xie(sigma*k, 1, 0, (nu-mu)/sigma, xi_guess_rescaled)*k*sigma + mu*k;
    count = count+1;
end
%% Figures
close all
figure
plot(kplot, imag(omega_spectral)); hold on
plot(kplot, imag(omega_xie),'.-');
plot(kplot, imag(omega_xie_rescaled));
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
title('IncompleteMaxwellian - $k$ vs. $\gamma$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Xie Rescale','location','Best')

figure
plot(kplot, real(omega_spectral)); hold on
plot(kplot, real(omega_xie),'.-');
plot(kplot, real(omega_xie_rescaled));
title('IncompleteMaxellian - $k$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Xie Rescale','Location','Best')