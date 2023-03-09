kplot = 0.1:0.005:1;
count = 1;

kappa = 1;
theta = [0.5, 0.5];
mu = [0.5, 3.5];
beta = 0.7;

moduli = zeros(1,length(kplot));

for k=kplot
    init_guess = Vlasov_1D_linearized_Steve_v3_Kappa(k, theta, mu, beta, kappa);

    zr = real(init_guess);
    zi = imag(init_guess);
    disp_relation = 1 + sqrt(2)*32/k*1i*theta(1)^3*beta*(  ...
    k*(sqrt(2)*1i*zi - 3*1i*k*theta(1) - sqrt(2)*k*mu(1) + sqrt(2)*zr)/...
    (8*theta(1)^3*(-2*zi + sqrt(2)*k*theta(1) - 2*1i*k*mu(1) + 2*1i*zr)^3) + ...
    k^4*(-1i*zi + k*mu(1) - zr)/...
    (2*zi^2 - k^2*(theta(1)^2+2*mu(1)^2) + 4*1i*zi*(k*mu(1)-zr) + 4*k*mu(1)*zr - 2*zr^2 )^3  )...
    + sqrt(2)*32/k*1i*theta(2)^3*(1-beta)*(  ...
    k*(sqrt(2)*1i*zi - 3*1i*k*theta(2) - sqrt(2)*k*mu(2) + sqrt(2)*zr)/...
    (8*theta(2)^3*(-2*zi + sqrt(2)*k*theta(2) - 2*1i*k*mu(2) + 2*1i*zr)^3) + ...
    k^4*(-1i*zi + k*mu(2) - zr)/...
    (2*zi^2 - k^2*(theta(2)^2+2*mu(2)^2) + 4*1i*zi*(k*mu(2)-zr) + 4*k*mu(2)*zr - 2*zr^2 )^3  );

    moduli(count) = abs(disp_relation-init_guess);

end

figure
plot(kplot,moduli);