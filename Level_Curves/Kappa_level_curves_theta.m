thetaplot = 0.25:0.01:1.75;
count=1;
initial_guesses = zeros(1,length(kplot));
gamma_xie = zeros(1,length(kplot));
gamma_root = zeros(1,length(kplot));

kappa = 1;
k = 0.5;
mu = 0;

tic;

for theta=thetaplot
    init_guess = Vlasov_1D_linearized_Steve_v3_Kappa(k, theta, mu, kappa);
    initial_guesses(count) = init_guess;
    gamma_xie(count) = Kappa_Disp_Using_Xie(k, theta, mu, kappa, init_guess);
    gamma_root(count) = dielectric(k,theta,theta,mu,mu,1,init_guess);
    count = count+1;
end

toc

% Create Plots
figure
plot(thetaplot, imag(initial_guesses))
hold on
plot(thetaplot, imag(gamma_xie))
plot(thetaplot, imag(gamma_root))
xlabel('$\theta$','Interpreter','latex','FontSize',16)
ylabel('$\gamma(\theta)$','Interpreter','latex','FontSize',16)
title(['Kappa Distribution ($\kappa=',num2str(kappa),'$)'],'Interpreter','latex')
legend('Spectral Method','Dispersion Relation with Xie Code','Root Finding with dielectric.m')

figure
plot(thetaplot, real(initial_guesses));
hold on
plot(thetaplot, real(gamma_xie));
plot(thetaplot, real(gamma_root));
xlabel('$\theta$','Interpreter','latex','FontSize',16)
ylabel('$\Omega(\theta)$','Interpreter','latex','FontSize',16)
title('Real Part of Initial Guess and Xie')
legend('Spectral Method','Dispersion Relation with Xie Code','Root Finding with dielectric.m')
