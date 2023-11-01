clear; clc; close all
nu0 = -2; nuf = 2;
nuplot = nu0:(nuf-nu0)/99:nuf;
count = 1;
omega_spectral = zeros(1,length(nuplot));
omega_xie = zeros(1,length(nuplot));
omega_xie_rescaled = zeros(1,length(nuplot));

sigma = 1;
mu = 0;
M = 10;
k = 0.5;
% vplot = -2.5:0.005:2.5;
% f0 = (atan(M*(vplot-nu))/pi+1/2).*exp(-(vplot-mu).^2/sigma^2)/sqrt(pi*sigma^2);
% plot(vplot,f0)

for nu=nuplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0, nu, M);
    xi_guess = (init_guess+mu*k)/k;
    xi_guess_rescaled = init_guess/(sigma*k);
    
    %If constant initial guess is plugged in Xie finds well continuous
    %looking solutions. no major jumps
    %xi_guess = (1.2-0.1i)/k;
    %xi_guess_rescaled = ((1.2-0.1i)-mu*k)/(sigma*k);

    omega_spectral(count) = init_guess + mu*k;
    omega_xie(count) = IncompleteMax_Disp_Using_Xie(k, sigma, mu, nu, xi_guess)*k;
    omega_xie_rescaled(count)= IncompleteMax_Disp_Using_Xie(sigma*k, 1, 0, (nu-mu)/sigma, xi_guess_rescaled)*k*sigma + mu*k;
    count = count+1;
end
%% Figures

% Imaginary part
figure
plot(nuplot, imag(omega_spectral)); hold on
plot(nuplot, imag(omega_xie),'.-');
plot(nuplot, imag(omega_xie_rescaled));
xlabel('$\nu$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\nu)$','Interpreter','latex','FontSize',16)
title('Incomplete Maxwellian - $\nu$ vs. $\gamma$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Xie Rescale','location','Best')

% Real plot
figure
plot(nuplot, real(omega_spectral)); hold on
plot(nuplot, real(omega_xie),'.-');
plot(nuplot, real(omega_xie_rescaled));
title('Incomplete Maxwellian - $\nu$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$\nu$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\nu)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Xie Rescale','location','Best')