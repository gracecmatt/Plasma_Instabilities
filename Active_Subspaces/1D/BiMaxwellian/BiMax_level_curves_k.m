clear; clc; 
k0 = 0.2; kf = 0.8;
kplot = k0:(kf-k0)/199:kf;
% kplot = kf:-(kf-k0)/199:k0;
count = 1;
initial_guesses = zeros(1,length(kplot));
omega_xie = zeros(1,length(kplot));
omega_zeta = zeros(1,length(kplot));
omega_xie_rescaled = zeros(1,length(kplot));

sigma1 = 1;
sigma2 = 1;
mu1 = 0;
mu2 = mu1+2;
beta = 0.9;

for k=kplot
    % init_guess = dispersion_growthrate_BiMax([k,sigma1,sigma2,mu1,mu2,beta],1-0.01*1i)*k - mu1*k;
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma1, sigma2, 0, mu2-mu1, beta); %\tilde{Omega}+igamma
    % init_guess = analyticAPPROX(k, sigma1, sigma2, 0, mu2-mu1, beta);
    initial_guesses(count) = init_guess + mu1*k; %Omega+igamma
    
    xi = (init_guess + mu1*k)/k; %xi=omega/k
    xi_scaled = (init_guess)/(sigma1*k);

    omega_xie(count) = BiMaxwellian_Disp_Using_Xie(k, sigma1, sigma2, mu1, mu2, beta, xi)*k; %omega=xi*k
    omega_xie_rescaled(count) = BiMaxwellian_Disp_Using_Xie(k*sigma1, 1, sigma2/sigma1, 0, (mu2-mu1)/sigma1, beta, xi_scaled)*k*sigma1 + mu1*k; % omega = xi*k*sigma1 + mu1*k
    omega_zeta(count) = dispersion_growthrate_BiMax([k,sigma1,sigma2,mu1,mu2,beta],xi)*k;
    count = count+1;
end

%% Figures
close all
v=linspace(mu1-5,mu1+5,1000);
f_gauss = @(x) exp(-x.^2)/sqrt(pi);
f0=beta*f_gauss((v-mu1)/sigma1)/sigma1 + (1-beta)*f_gauss((v-mu2)/sigma2)/sigma2;
figure; plot(v,f0,'linewidth',2); title('Plot of Velocity Distribution $f_0(v)$','Interpreter','latex','FontSize',16)

txt1 = ['$\mu_1$ = ',num2str(mu1),', $\mu_2$ = ', num2str(mu2)];
txt2 = ['$\sigma_1$ = ',num2str(sigma1),', $\sigma_2$ = ', num2str(sigma2)];
txt3 = ['$\beta$ = ',num2str(beta)];
txt = {txt1,txt2,txt3};

figure
plot(kplot, imag(initial_guesses),'.-'); hold on
plot(kplot, imag(omega_xie),'.-');
plot(kplot, imag(omega_xie_rescaled));
plot(kplot, imag(omega_zeta));
title('Bi-Maxwellian - $k$ vs. $\gamma$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method','Xie Code','Rescaled Xie','Zetaf Disp.','location','Best')
xL=xlim; yL=ylim;
text(kplot(floor(end/2)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

figure
plot(kplot, real(initial_guesses),'.-'); hold on
plot(kplot, real(omega_xie),'.-');
plot(kplot, real(omega_xie_rescaled));
plot(kplot, real(omega_zeta));
title('Bi-Maxwellian - $k$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method','Xie Code','Rescaled Xie','Zetaf Disp.','Location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(kplot(2)-kplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)