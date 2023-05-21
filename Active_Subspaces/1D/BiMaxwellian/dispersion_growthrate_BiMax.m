function gamma = dispersion_growthrate_BiMax(p)
% Calculating growth rate (imaginary component of angular frequency) using
% explicit root finding method for dispersion relation.
% ---Requires "zetaf.m" - representation of Z function
% ---Input p is an 8-parameter array in the form p=[alpha k mu sigma]

% BiMaxwellian/Bump-on-Tail
% %
% beta_1 = p(1) ; beta_2 =  1-beta_1;
% sigma_1 = p(2) ; sigma_2 = p(4) ;
% mu_1 =  p(3); mu_2  = p(5) ;
% 
% k = p(6);

k = p(1);

mu_1 =  p(2); mu_2  = p(3) ;
sigma_1 = p(4) ; sigma_2 = p(5) ;
beta_1 = p(6) ; beta_2 =  1-beta_1;

guess_real = 1;%1; %p(8)
guess_imag = 0.5;%0.5; %p(9)


A1 = @(zr, zi) (1/(sqrt(2*sigma_1)))*(zr + 1i*zi- mu_1);
A2 = @(zr, zi) (1/(sqrt(2*sigma_2)))*(zr + 1i*zi- mu_2);

F = @(z) [real(1+(beta_1/(sigma_1*k^2))*(1+A1(z(1),z(2))*zetaf(A1(z(1),z(2))))+(beta_2/(sigma_2*k^2))*(1+A2(z(1),z(2))*zetaf(A2(z(1),z(2))))); 
            imag(1+(beta_1/(sigma_1*k^2))*(1+A1(z(1),z(2))*zetaf(A1(z(1),z(2))))+(beta_2/(sigma_2*k^2))*(1+A2(z(1),z(2))*zetaf(A2(z(1),z(2)))))];

options = optimset('Display','off');
z_r1 = fsolve(F, [guess_real, guess_imag], options);
%gamma = k*(z_r1(1) +1i*z_r1(2)); % complex-valued angular frequency
gamma = k*z_r1(2); %imaginary component of angular frequency for GSA

       
