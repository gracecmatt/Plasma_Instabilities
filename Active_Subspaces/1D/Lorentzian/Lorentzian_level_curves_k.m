%% Lorentzian
clear; clc; 
k0 = 0.5; kf = 1.5;
kplot = k0:(kf-k0)/499:kf;
count = 1;
initial_guesses = zeros(1,length(kplot));
omega_xie = zeros(1,length(kplot));
% omega_xie_rescaled = zeros(1,length(kplot));

sigma = 1;
mu = 900;
exactReal = mu*kplot+1; % other solution: mu*kplot-1
exactImag = -sigma*kplot;

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma, 0); %tilde{Omega}+igamma
    initial_guesses(count) = init_guess+mu*k; %Omega+igamma

    xi = (init_guess+mu*k)/k; 
    % xi_rescaled = init_guess/(sigma*k);

    omega_xie(count) = Lorentzian_Disp_Using_Xie(k, sigma, mu, xi)*k; %omega=xi*k
    % omega_xie_rescaled(count) = Lorentzian_Disp_Using_Xie(k*sigma, 1, 0, xi_rescaled)*sigma*k + mu*k; %omega=xi*sigma*k+mu*k
    count = count+1;
end

%% Figures
close all
txt = ['$\mu$ = ',num2str(mu),', $\sigma$ = ', num2str(sigma)];

figure
plot(kplot, imag(initial_guesses),'.-'); hold on
plot(kplot, imag(omega_xie),'.-');
plot(kplot, imag(omega_xie_rescaled),'.-');
plot(kplot, exactImag,'k');
title('Lorentzian - $k$ vs. $\gamma$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','Exact Solution','location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(kplot(2)-kplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

figure
plot(kplot, real(initial_guesses),'.-'); hold on
plot(kplot, real(omega_xie),'.-');
plot(kplot, real(omega_xie_rescaled),'.-');
plot(kplot, exactReal, 'k')
title('Lorentzian - $k$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','Exact Solution','Location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(kplot(2)-kplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

%% Error Analysis
% L2 error = sum( (y_exact - y_sample).^2 )

% real part: Omega = mu*k+1
L2err.spectral(1) = sum( (exactReal-(real(initial_guesses))).^2 );
L2err.xie(1) = sum( (exactReal-(real(omega_xie))).^2 );
L2err.xie_rescaled(1) = sum( (exactReal-(real(omega_xie_rescaled))).^2 );

% imaginary part: gamma = -sigma*k
L2err.spectral(2) = sum( (exactImag-(imag(initial_guesses))).^2 );
L2err.xie(2) = sum( (exactImag-(imag(omega_xie))).^2 );
L2err.xie_rescaled(2) = sum( (exactImag-(imag(omega_xie_rescaled))).^2 );

Error = struct2table(L2err)