function gamma = Kappa_Disp_Relation(k, theta, mu, init_guess)

% Kappa = 1 throughout

% Real and imaginary parts
guess_real = real(init_guess);
guess_imag = imag(init_guess);

% A = @(zr, zi) zr + zi*1i + 1/sqrt(2)*k*theta*1i; %fixed GM 12/9/22
% F = @(z) [real( 1 - (sqrt(2)*k*theta*1i + A(z(1),z(2))) / (A(z(1),z(2)))^3 )
%           imag( 1 - (sqrt(2)*k*theta*1i + A(z(1),z(2))) / (A(z(1),z(2)))^3 )];

% With mu; added by GM 12/23/22
A = @(zr, zi) 1 + sqrt(2)*32/k*1i*theta^3*( ...
k*(sqrt(2)*1i*zi - 3*1i*k*theta - sqrt(2)*k*mu + sqrt(2)*zr)/...
(8*theta^3*(-2*zi + sqrt(2)*k*theta - 2*1i*k*mu + 2*1i*zr)^3) + ...
k^4*(-1i*zi + k*mu - zr)/...
(2*zi^2 - k^2*(theta^2+2*mu^2) + 4*1i*zi*(k*mu-zr) + 4*k*mu*zr - 2*zr^2 )^3 );
F = @(z) [real( A(z(1),z(2)) )
          imag( A(z(1),z(2)) )];

options = optimoptions('fsolve','Display','off');
z_r1 = fsolve(F, [guess_real, guess_imag], options);
gamma = z_r1(2); %imaginary component of angular frequency for GSA

end

