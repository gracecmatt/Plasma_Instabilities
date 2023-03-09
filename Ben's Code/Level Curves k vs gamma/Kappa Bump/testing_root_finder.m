k = 0.33;
mu1 = 0;
mu2 = 4;
theta1 = 1;
theta2 = 1;
beta = 0.8;

omega_guess = Vlasov_1D_linearized_Steve_v3_Kappa(k, [theta1, theta2], [mu1, mu2], beta, kappa);
omega_real = real(omega_guess);
gamma = imag(omega_guess);

disp(['Expected Value: ', gamma])

[x, y] = meshgrid(omega_real-0.2:0.01:omega_real+0.2, gamma-0.2:0.01:gamma+0.2);
complex_grid = x + 1i*y;
roots = zeros(size(complex_grid));

[rows, cols] = size(complex_grid);
for i=1:rows
    for j=1:cols
        init_guess = complex_grid(i,j);
        gamma = dielectric(k, theta1, theta2, mu1, mu2, theta, init_guess);
        roots(i,j) = gamma;
    end
end

figure
xlabel('Real')
ylabel('Imag')
title('Roots')
surf(x, y, roots);
figure
xlabel('Real')
ylabel('Imag')
title('Difference of initial guess and root')
surf(x, y, abs(roots-gamma))