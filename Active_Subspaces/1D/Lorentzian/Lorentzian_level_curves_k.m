%% Lorentzian
clear; clc; 
kplot = linspace(0.25,0.75);
count = 1;
spectral_guess = zeros(1,length(kplot));
omega = zeros(1,length(kplot));
omega_rescaled = zeros(1,length(kplot));

sigma = 1;
mu = 100;

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0); %tilde{Omega}+igamma
    
    xi_guess = (init_guess+mu*k)/k; %xi=omega/k
    xi_guess_rescaled = init_guess/(sigma*k);

    spectral_guess(count) = init_guess+mu*k; %Omega+igamma
    omega(count) = Lorentzian_Disp_Using_Xie(k, sigma, mu, xi_guess)*k; %omega=xi*k
    omega_rescaled(count) = Lorentzian_Disp_Using_Xie(k*sigma, 1, 0, xi_guess_rescaled)*sigma*k + mu*k; %omega=xi*sigma*k+mu*k
    
    count = count+1;
end

% exact solution
omega_exact = mu*kplot+1 + 1i.*(-sigma*kplot);

%% Figures
close all
txt = ['$\mu$ = ',num2str(mu),', $\sigma$ = ', num2str(sigma)];

figure
plot(kplot, imag(spectral_guess),'.-'); hold on
plot(kplot, imag(omega),'.-');
plot(kplot, imag(omega_rescaled),'.-');
plot(kplot, imag(omega_exact),'k');
title('Lorentzian $\gamma(k)$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','Exact Solution','location','East')
xL=xlim; yL=ylim;
text(median(kplot),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

figure
plot(kplot, real(spectral_guess),'.-'); hold on
plot(kplot, real(omega),'.-');
plot(kplot, real(omega_rescaled),'.-');
plot(kplot, real(omega_exact), 'k')
title('Lorentzian $\Omega(k)$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','Exact Solution','Location','South')
xL=xlim; yL=ylim;
text(xL(1)+(kplot(2)-kplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

%% Error Analysis
% L2 error = sqrt( sum( (y_exact - y_sample).^2 ) )
% results are given as two component vectors of [real l2 error, imag l2 error]

% real part
L2err.spectral(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(spectral_guess)) ).^2 ));
L2err.xie(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(omega)) ).^2 ));
L2err.xie_rescaled(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(omega_rescaled)) ).^2 ));

% imaginary part
L2err.spectral(2) = sqrt(sum( (imag(omega_exact)-(imag(spectral_guess))).^2 ));
L2err.xie(2) = sqrt(sum( (imag(omega_exact)-(imag(omega))).^2 ));
L2err.xie_rescaled(2) = sqrt(sum( (imag(omega_exact)-(imag(omega_rescaled))).^2 ));

Error = struct2table(L2err)