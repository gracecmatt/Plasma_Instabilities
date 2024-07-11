function xi = Max_dielectric(p,init_guess)
% Calculating growth rate (imaginary component of angular frequency) using
% explicit root finding method for dispersion relation.
% ---Requires "zetaf.m" - representation of Z function
% % ---Input p is an 8-parameter array in the form p=[alpha k mu sigma]

% NOTE: check F because I changed sigma in it so that it is consistent with
% the function we are using elsewhere.

% BiMaxwellian/Bump-on-Tail
k = p(1);
sigma = p(2);
mu = p(3);

guess_real = real(init_guess); 
guess_imag = imag(init_guess);

A = @(zr, zi) 1/sigma*(zr + 1i*zi- mu);

F = @(z) [real(1+(2/(sigma^2*k^2))*(1+A(z(1),z(2))*zetaf(A(z(1),z(2))))); 
            imag(1+(2/(sigma^2*k^2))*(1+A(z(1),z(2))*zetaf(A(z(1),z(2)))))];

options = optimset('Display','off');
z_r1 = fsolve(F, [guess_real, guess_imag], options);
xi = abs(z_r1(1)) + 1i*z_r1(2); % complex-valued angular frequency ==(omega/k)