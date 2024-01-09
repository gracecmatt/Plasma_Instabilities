function xi = BiMax_dielectric(p,init_guess)
% Calculating growth rate (imaginary component of angular frequency) using
% explicit root finding method for dispersion relation.
% ---Requires "zetaf.m" - representation of Z function
% ---Input p is an 8-parameter array in the form p=[alpha k mu sigma]

% NOTE: check F because I changed sigma in it so that it is consistent with
% the function we are using elsewhere.

% BiMaxwellian/Bump-on-Tail
k = p(1);
sigma_1 = p(2);         sigma_2 = p(3);
mu_1 =  p(4);           mu_2  = p(5);
beta_1 = p(6);          beta_2 =  1-beta_1;

guess_real = real(init_guess); 
guess_imag = imag(init_guess);

A1 = @(zr, zi) 1/sigma_1*(zr + 1i*zi- mu_1);
A2 = @(zr, zi) 1/sigma_2*(zr + 1i*zi- mu_2);

F = @(z) [real(1+(2*beta_1/(sigma_1^2*k^2))*(1+A1(z(1),z(2))*zetaf(A1(z(1),z(2))))+(2*beta_2/(sigma_2^2*k^2))*(1+A2(z(1),z(2))*zetaf(A2(z(1),z(2))))); 
            imag(1+(2*beta_1/(sigma_1^2*k^2))*(1+A1(z(1),z(2))*zetaf(A1(z(1),z(2))))+(2*beta_2/(sigma_2^2*k^2))*(1+A2(z(1),z(2))*zetaf(A2(z(1),z(2)))))];

options = optimset('Display','off');
z_r1 = fsolve(F, [guess_real, guess_imag], options);
xi = abs(z_r1(1)) + 1i*z_r1(2); % complex-valued angular frequency ==(omega/k)