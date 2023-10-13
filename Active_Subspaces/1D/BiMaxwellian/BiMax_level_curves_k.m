clear; clc;
kplot = linspace(0.25,0.75);
count = 1;
spectral_guess = zeros(1,length(kplot));
omega = zeros(1,length(kplot));
omega_rescaled = zeros(1,length(kplot));
omega_dielectric = zeros(1,length(kplot));

sigma1 = 5;
sigma2 = 1;
mu1 = 0;
mu2 = mu1+4;
beta = 0.9;

% eyeball-check the distribution
v=linspace(mu1-10,mu1+10,1000);
f_gauss = @(x) exp(-x.^2)/sqrt(pi);
f0= @(v) beta*f_gauss((v-mu1)/sigma1)/sigma1 + (1-beta)*f_gauss((v-mu2)/sigma2)/sigma2;
figure; plot(v,f0(v),'linewidth',2); title('Plot of Velocity Distribution $f_0(v)$','Interpreter','latex','FontSize',16)
% check normalization
integral(f0,-Inf,Inf); pause(0.1);

tic;
for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v4(k, sigma1, sigma2, 0, mu2-mu1, beta); %\tilde{Omega}+igamma

    xi_guess = (init_guess + mu1*k)/k; %xi=omega/k
    xi_guess_rescaled = (init_guess)/(sigma1*k); 

    spectral_guess(count) = init_guess + mu1*k; %omega=Omega+igamma
    omega(count) = BiMaxwellian_Disp_Using_Xie(k, sigma1, sigma2, mu1, mu2, beta, xi_guess)*k; %omega=xi*k
    omega_rescaled(count) = BiMaxwellian_Disp_Using_Xie(k*sigma1, 1, sigma2/sigma1, 0, (mu2-mu1)/sigma1, beta, xi_guess_rescaled)*k*sigma1 + mu1*k; % omega = xi*k*sigma1 + mu1*k
    omega_dielectric(count) = BiMax_dispersion([k,sigma1,sigma2,mu1,mu2,beta],xi_guess)*k; 
    
    count = count+1;
end
toc

%% Figures
close all
% % plot distribution
% v=linspace(mu1-10,mu1+10,1000);
% f_gauss = @(x) exp(-x.^2)/sqrt(pi);
% f0=@(v)beta*f_gauss((v-mu1)/sigma1)/sigma1 + (1-beta)*f_gauss((v-mu2)/sigma2)/sigma2;
% figure; plot(v,f0(v),'linewidth',2); title('Plot of Velocity Distribution $f_0(v)$','Interpreter','latex','FontSize',16)

txt1 = ['$\mu_1$ = ',num2str(mu1),', $\mu_2$ = ', num2str(mu2)];
txt2 = ['$\sigma_1$ = ',num2str(sigma1),', $\sigma_2$ = ', num2str(sigma2)];
txt3 = ['$\beta$ = ',num2str(beta)];
txt = {txt1,txt2,txt3};

figure
plot(kplot, imag(spectral_guess),'.-'); hold on
plot(kplot, imag(omega),'.-');
plot(kplot, imag(omega_rescaled),'.-');
plot(kplot, imag(omega_dielectric),'o-');
title('BiMaxwellian $\gamma(k)$ ','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method','Xie Code','Rescaled Xie','Zetaf Disp.','location','East')
xL=xlim; yL=ylim;
text(median(kplot),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

figure
plot(kplot, real(spectral_guess),'.-'); hold on
plot(kplot, real(omega),'.-');
plot(kplot, real(omega_rescaled),'.-');
plot(kplot, real(omega_dielectric),'o-');
title('BiMaxwellian $\Omega(k)$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method','Xie Code','Rescaled Xie','Zetaf Disp.','Location','East')
xL=xlim; yL=ylim;
text(xL(1)+(kplot(2)-kplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

%% Error Analysis (using dielectric function as "true" solution)
% L2 error = sqrt( sum( (y_exact - y_sample).^2 ) )
% results are given as two component vectors of [real l2 error, imag l2 error]

% real part
L2err.spectral(1) = sqrt(sum( ( abs(real(omega_dielectric))-abs(real(spectral_guess)) ).^2 ));
L2err.xie(1) = sqrt(sum( ( abs(real(omega_dielectric))-abs(real(omega)) ).^2 ));
L2err.xie_rescaled(1) = sqrt(sum( ( abs(real(omega_dielectric))-abs(real(omega_rescaled)) ).^2 ));

% imaginary part
L2err.spectral(2) = sqrt(sum( (imag(omega_dielectric)-(imag(spectral_guess))).^2 ));
L2err.xie(2) = sqrt(sum( (imag(omega_dielectric)-(imag(omega))).^2 ));
L2err.xie_rescaled(2) = sqrt(sum( (imag(omega_dielectric)-(imag(omega_rescaled))).^2 ));

Error = struct2table(L2err)