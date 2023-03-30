function gamma = Kappa_Disp_Relation_BumponTail(k, theta, mu, beta, init_guess)

% Kappa = 1 throughout
% mu and theta are 2x1 vectors with mu(2),theta(2) = mean velocity/std of bump

% Real and imaginary parts
guess_real = real(init_guess);
guess_imag = imag(init_guess);

% With mu and bump on tail; added by GM 12/23/22
A = @(zr, zi) 1 + sqrt(2)*32/k*1i*theta(1)^3*beta*(  ...
k*(sqrt(2)*1i*zi - 3*1i*k*theta(1) - sqrt(2)*k*mu(1) + sqrt(2)*zr)/...
(8*theta(1)^3*(-2*zi + sqrt(2)*k*theta(1) - 2*1i*k*mu(1) + 2*1i*zr)^3) + ...
k^4*(-1i*zi + k*mu(1) - zr)/...
(2*zi^2 - k^2*(theta(1)^2+2*mu(1)^2) + 4*1i*zi*(k*mu(1)-zr) + 4*k*mu(1)*zr - 2*zr^2 )^3  )...
+ sqrt(2)*32/k*1i*theta(2)^3*(1-beta)*(  ...
k*(sqrt(2)*1i*zi - 3*1i*k*theta(2) - sqrt(2)*k*mu(2) + sqrt(2)*zr)/...
(8*theta(2)^3*(-2*zi + sqrt(2)*k*theta(2) - 2*1i*k*mu(2) + 2*1i*zr)^3) + ...
k^4*(-1i*zi + k*mu(2) - zr)/...
(2*zi^2 - k^2*(theta(2)^2+2*mu(2)^2) + 4*1i*zi*(k*mu(2)-zr) + 4*k*mu(2)*zr - 2*zr^2 )^3  );

F = @(z) [real( A(z(1),z(2)) )
          imag( A(z(1),z(2)) )];

options = optimoptions('fsolve','Display','none');
z_r1 = fsolve(F, [guess_real, guess_imag], options);

if(abs(z_r1(2)-guess_imag) > 0.02)
    gamma = guess_imag;
else
    gamma = z_r1(2);
end

end

