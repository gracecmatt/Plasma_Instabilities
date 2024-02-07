%% Non-Integer Kappa Comparison - Mace & Hellburg (1995) 
clear; clc;

kappa = 1.63;
% assuming: mu = 0; sigma = 1;
z = linspace(0.1,1.1,24)+1i*linspace(0.1,1.1,24); % sample random complex z values


% ========= compute Z, Zp using Fourier series approximation ==========
F = ['(gamma(',num2str(kappa,16),')/gamma(',num2str(kappa,16),'-0.5)/sqrt(pi*',...
    num2str(kappa,16),'))*(1+v.^2/',num2str(kappa,16),').^(-',num2str(kappa,16),'-1)'];
Fn = 0; % 0 = option to define your own F
Nfourier = 400; % number of Fourier coefficients to take
[Zp,Z] = zetaph(z, Fn, F, Nfourier); % dielectric function


% ========= compute Z, Zp using Gauss hypergeometric function ==========
x = 1/2*(1-z/(1i*sqrt(kappa)));
Zgauss = 1i*(kappa+1/2)*(kappa-1/2)/(kappa^(3/2)*(kappa+1))*hypergeom( [1, 2*kappa+2], kappa+2, x );
Zpgauss = -(kappa+1/2)*(kappa-1/2)/(kappa^2*(kappa+2))*hypergeom( [2, 2*kappa+3], kappa+3, x );


% ========= difference between algorithms ==========
max(abs(Z-Zgauss))
max(abs(Zp-Zpgauss))

