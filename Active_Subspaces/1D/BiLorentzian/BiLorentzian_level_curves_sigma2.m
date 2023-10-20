%% BiLorentzian
clear; clc;
s2plot = linspace(0.3,0.9,100);
count = 1;
spectral_guess = zeros(1,length(s2plot));
omega = zeros(1,length(s2plot));
omega_rescaled = zeros(1,length(s2plot));
omega_exact = zeros(1,length(s2plot));

k = 0.5;
sigma2 = 0.6;
sigma1 = 1;
mu1 = 50;
mu2 = mu1+6;
beta = 0.9;

% % view the distribution you are about to study
v = linspace(mu1-10,mu1+10,2000);
f0 = @(v) beta/pi*sigma1./((v-mu1).^2+sigma1^2) + (1-beta)/pi*sigma2./((v-mu2).^2+sigma2^2);
figure; plot(v,f0(v),'linewidth',2); title('Plot of Velocity Distribution $f_{BL}(v)$','Interpreter','latex','FontSize',16)
% check normalization
integral(f0,-Inf,Inf); pause(0.1);

for sigma2=s2plot
      % returns abs(tilde{Omega}) + igamma
    % init_guess = analyticAPPROX(k,sigma1,sigma2,0,mu2-mu1,beta);
    init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma1,sigma2,0,mu2-mu1,beta); %tilde{Omega}+igamma

    xi_guess = (init_guess+mu1*k)/k;
    xi_guess_rescaled = init_guess/(sigma1*k);

    spectral_guess(count) = init_guess+mu1*k; %Omega+igamma
      % returns tilde{Omega} + igamma
    omega(count) = BiLorentzian_Disp_Using_Xie(k,sigma1,sigma2,mu1,mu2,beta,xi_guess)*k; %omega=xi*k
    omega_rescaled(count) = BiLorentzian_Disp_Using_Xie(k*sigma1,1,sigma2/sigma1,0,(mu2-mu1)/sigma1,beta,xi_guess_rescaled)*sigma1*k + mu1*k; %omega=xi*sigma1*k+mu1*k
      % returns abs(tilde{Omega})+mu1*k + igamma
    omega_exact(count) = BiLorentzian_Solution(k,sigma1,sigma2,mu1,mu2,beta);

    count = count+1;
end

%% Figures
close all
txt1 = ['$\mu_1$ = ',num2str(mu1),', $\mu_2$ = ', num2str(mu2)];
txt2 = ['$\sigma_1$ = ',num2str(sigma1),', $\sigma_2$ = ', num2str(sigma2)];
txt3 = ['$\beta$ = ',num2str(beta)];
txt = {txt1,txt2,txt3};

figure
plot(s2plot, imag(spectral_guess),'.-'); hold on
plot(s2plot, imag(omega),'.-');
plot(s2plot, imag(omega_rescaled),'.-');
plot(s2plot, imag(omega_exact),'k');
title('BiLorentzian $\gamma(\sigma_2)$','Interpreter','latex','FontSize',16)
xlabel('$\sigma_2$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\sigma_2)$','Interpreter','latex','FontSize',16)
legend('Spectral Method','Xie', 'Rescaled Xie','Exact Solution','location','East')
xL=xlim; yL=ylim;
text(median(s2plot),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

figure
plot(s2plot, abs(real(spectral_guess)),'.-'); hold on
plot(s2plot, abs(real(omega)),'.-');
plot(s2plot, abs(real(omega_rescaled)),'.-');
plot(s2plot, abs(real(omega_exact)), 'k');
title('BiLorentzian $\Omega(\sigma_2)$','Interpreter','latex','FontSize',16)
xlabel('$\sigma_2$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\sigma_2)$','Interpreter','latex','FontSize',16)
legend('Spectral Method','Xie','Rescaled Xie','Exact Solution','Location','East')
xL=xlim; yL=ylim;
text(median(s2plot),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

%% Error Analysis
% L2 error = sqrt( sum( (y_exact - y_sample).^2 ) )
% results are given as two component vectors of [real l2 error, imag l2 error]

% NOTE: if error for real part of xie_rescaled is high, adjust spectral
% code to be more accurate and error should decrease.

% real part
L2err.spectral(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(spectral_guess)) ).^2 ));
L2err.xie(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(omega)) ).^2 ));
L2err.xie_rescaled(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(omega_rescaled)) ).^2 ));

% imaginary part
L2err.spectral(2) = sqrt(sum( (imag(omega_exact)-(imag(spectral_guess))).^2 ));
L2err.xie(2) = sqrt(sum( (imag(omega_exact)-(imag(omega))).^2 ));
L2err.xie_rescaled(2) = sqrt(sum( (imag(omega_exact)-(imag(omega_rescaled))).^2 ));

Error = struct2table(L2err)