clear; clc;
kplot = linspace(0.25,0.75);
count = 1;
spectral_guess = zeros(1,length(kplot));
omega = zeros(1,length(kplot));
omega_rescaled = zeros(1,length(kplot));
omega_dielectric = zeros(1,length(kplot));

kappa = 2;
sigma = 1;
mu = 10;

% eyeball-check the distribution
v=linspace(mu-10,mu+10,1000);
C1=(pi*sigma^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
f0 = @(v) C1*(1+(v-mu).^2/((kappa-1.5)*sigma^2)).^(-kappa);
figure; plot(v,f0(v),'linewidth',2); title('Plot of Velocity Distribution $f_0(v)$','Interpreter','latex','FontSize',16)
% check normalization
integral(f0,-Inf,Inf)

tic;
for k=kplot

    init_guess = Vlasov_1D_linearized_Steve_v4_Kappa(k, sigma, mu, kappa); %tilde{omega}=tilde{Omega}+igamma
    
    omega_guess = init_guess+mu*k; %omega=Omega+igamma
    xi_guess = (init_guess+mu*k)/k; %xi=omega/k
    xi_guess_rescaled = init_guess/(sigma*k); 

    spectral_guess(count) = init_guess + mu*k; %omega=Omega+igamma
    omega(count) = Kappa_Disp_Using_Xie(k,sigma,mu,kappa,xi_guess)*k; %omega=xi*k
    omega_rescaled(count) = Kappa_Disp_Using_Xie(k*sigma,1,0,kappa,xi_guess_rescaled)*sigma*k + mu*k; %omega=xi*sigma1*k+mu1*k
    omega_dielectric(count) = dielectric_kappa(k,sigma,0.5,mu,0,1,kappa,omega_guess); % This is not the correct function

    count = count+1;
end
toc

%% Figures 
close all;
% % plot distribution
% v=linspace(mu1-10,mu1+10,1000);
% C1=(pi*sigma1^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
% C2=(pi*sigma2^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
% f0= @(v) beta*C1*(1+(v-mu1).^2/((kappa-1.5)*sigma1^2)).^(-kappa) + ...
%    (1-beta)*C2*(1+(v-mu2).^2/((kappa-1.5)*sigma2^2)).^(-kappa);
% figure; plot(v,f0(v),'linewidth',2); title('Plot of Velocity Distribution $f_0(v)$','Interpreter','latex','FontSize',16)

txt1 = ['$\mu_1$ = ',num2str(mu)];
txt2 = ['$\sigma_1$ = ',num2str(sigma)];
txt3 = ['$\kappa$ = ',num2str(kappa)];
txt = {txt1,txt2,txt3};

figure
plot(kplot, imag(spectral_guess),'.-'); hold on
plot(kplot, imag(omega),'.-');
plot(kplot, imag(omega_rescaled),'.-');
% plot(kplot, imag(omega_dielectric),'o-');
    titletxt = ['BiKappa $\gamma(k)$, $\kappa=$',num2str(kappa)];
    title(titletxt,'Interpreter','latex','FontSize',16)
    xlabel('$k$','Interpreter','latex','FontSize',16)
    ylabel('$\gamma(k)$','Interpreter','latex','FontSize',16)
    legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','Dielectric Func.',Location='South')
xL=xlim; yL=ylim;
text(median(kplot),yL(2),txt,'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','latex','FontSize',12)

figure
plot(kplot, real(spectral_guess),'.-'); hold on
plot(kplot, real(omega),'.-');
plot(kplot, real(omega_rescaled),'.-');
plot(kplot, real(omega_dielectric),'o-');
    titletxt = ['BiKappa $\Omega(k)$, $\kappa=$',num2str(kappa)];
    title(titletxt,'Interpreter','latex','FontSize',16)
    xlabel('$k$','Interpreter','latex','FontSize',16)
    ylabel('$\Omega(k)$','Interpreter','latex','FontSize',16)
    legend('Spectral Method', 'Xie Root Finding','Rescaled Xie','Dielectric Func.',Location='East')
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