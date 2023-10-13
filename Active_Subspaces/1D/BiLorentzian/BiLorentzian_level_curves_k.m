%% Lorentzian
clear; clc; 
k0 = 0.5; % center k of interest
var = 0.50; % x100 = % variation
N = 50; % +1 = number of points
ki = (1-var)*k0; kf = (1+var)*k0; 
kplot = [linspace(ki,k0,N/2),linspace(k0+(kf-k0)/(N/2+1),kf,N/2+1)];

count = 1;
initial_guesses = zeros(1,length(kplot));
omega_xie = zeros(1,length(kplot));
omega_xie_rescaled_1 = zeros(1,length(kplot));
omega_xie_rescaled_2 = zeros(1,length(kplot));
omega_exact = zeros(1,length(kplot));

sigma1 = 1;
sigma2 = 1;
mu1 = 0;
mu2 = mu1+6;
beta = 0.7;

% view the distribution you are about to study
v = linspace(mu1-8,mu1+8,2000);
f0 = @(v) beta/pi*sigma1./((v-mu1).^2+sigma1^2) + (1-beta)/pi*sigma2./((v-mu2).^2+sigma2^2);
figure; plot(v,f0(v),'linewidth',2); title('Plot of Velocity Distribution $f_0(v)$','Interpreter','latex','FontSize',16)
% check normalization
integral(f0,-Inf,Inf)

for k=kplot

    params(1) = k;
    params(2) = sigma1;
    params(3) = sigma2;
    params(4) = mu1;
    params(5) = mu2;
    params(6) = beta;

    % init_guess = Vlasov_1D_linearized_Steve_v4(k,sigma1,sigma2,0,mu2-mu1,beta); %tilde{Omega}+igamma
    init_guess = Vlasov_1D_linearized_Steve_v4(params(1),params(2),params(3),0,params(5)-params(4),params(6)); %tilde{Omega}+igamma
    % init_guess = BiLorentzian_Solution(k,sigma1,sigma2,mu1,mu2,beta) - mu1*k;
    initial_guesses(count) = init_guess+params(4)*params(1); %Omega+igamma

    init_guess_1 = abs(real(init_guess))+1i*imag(init_guess);
    init_guess_2 = -abs(real(init_guess))+1i*imag(init_guess);

    xi = (init_guess+mu1*params(1))/params(1); 
    xi_rescaled_1 = init_guess_1/(params(2)*params(1));
    xi_rescaled_2 = init_guess_2/(params(2)*params(1));

    omega_xie(count) = BiLorentzian_Disp_Using_Xie(k,sigma1,sigma2,mu1,mu2,beta,xi)*k; %omega=xi*k
    % omega_xie_rescaled(count) = BiLorentzian_Disp_Using_Xie(k*sigma1,1,sigma2/sigma1,0,(mu2-mu1)/sigma1,beta,xi_rescaled)*sigma1*k + mu1*k; %omega=xi*sigma*k+mu*k
    omega_xie_rescaled_1(count) = BiLorentzian_Disp_Using_Xie(params(1)*params(2),1,params(3)/params(2),0,(params(5)-params(4))/params(2),params(6),xi_rescaled_1)*params(2)*params(1) + params(4)*params(1); %omega=xi*sigma*k+mu*k
    omega_xie_rescaled_2(count) = BiLorentzian_Disp_Using_Xie(params(1)*params(2),1,params(3)/params(2),0,(params(5)-params(4))/params(2),params(6),xi_rescaled_2)*params(2)*params(1) + params(4)*params(1); %omega=xi*sigma*k+mu*k
  

    omega_exact(count) = BiLorentzian_Solution(params(1),params(2),params(3),params(4),params(5),params(6));
    count = count+1;
end

%% Figures
close all
txt1 = ['$\mu_1$ = ',num2str(mu1),', $\mu_2$ = ', num2str(mu2)];
txt2 = ['$\sigma_1$ = ',num2str(sigma1),', $\sigma_2$ = ', num2str(sigma2)];
txt3 = ['$\beta$ = ',num2str(beta)];
txt = {txt1,txt2,txt3};

figure
plot(kplot, imag(initial_guesses),'.-'); hold on
plot(kplot, imag(omega_xie),'.-');
plot(kplot, imag(omega_xie_rescaled_1),'o-');
plot(kplot, imag(omega_xie_rescaled_2),'.-');
plot(kplot, imag(omega_exact),'k');
% plot(kplot,ones(1,length(kplot))*imag(omega_xie_rescaled(floor(N/2))),'--k'); % reference guide for k=k0
title('Lorentzian - $k$ vs. $\gamma$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method','Xie', 'Rescaled Xie 1','Rescaled Xie 2','Exact Solution','location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(kplot(2)-kplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

figure
plot(kplot, real(initial_guesses),'.-'); hold on
plot(kplot, real(omega_xie),'.-');
plot(kplot, real(omega_xie_rescaled_1),'o-');
plot(kplot, real(omega_xie_rescaled_2),'.-');
plot(kplot, real(omega_exact), 'k');
title('Lorentzian - $k$ vs. $\Omega$','Interpreter','latex','FontSize',16)
xlabel('$k$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
legend('Spectral Method','Xie','Rescaled Xie 1','Rescaled Xie 2','Exact Solution','Location','Best')
xL=xlim; yL=ylim;
text(xL(1)+(kplot(2)-kplot(1)),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

%% Error Analysis
% L2 error = sqrt( sum( (y_exact - y_sample).^2 ) )
% results are given as two component vectors of [real l2 error, imag l2 error]

% NOTE: if error for real part of xie_rescaled is high, adjust spectral
% code to be more accurate and error should decrease.

% real part
L2err.spectral(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(initial_guesses)) ).^2 ));
L2err.xie(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(omega_xie)) ).^2 ));
L2err.xie_rescaled(1) = sqrt(sum( ( abs(real(omega_exact))-abs(real(omega_xie_rescaled)) ).^2 ));

% imaginary part
L2err.spectral(2) = sqrt(sum( (imag(omega_exact)-(imag(initial_guesses))).^2 ));
L2err.xie(2) = sqrt(sum( (imag(omega_exact)-(imag(omega_xie))).^2 ));
L2err.xie_rescaled(2) = sqrt(sum( (imag(omega_exact)-(imag(omega_xie_rescaled))).^2 ));

Error = struct2table(L2err)