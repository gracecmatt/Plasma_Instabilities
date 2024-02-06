function [omega] = AnalyticApproxKappa(k,sigma1,sigma2,v0,beta,kappa)
% plasma frequencies and total number of particles are normalized to one

    omegaRE = sqrt(1+3*k^2*sigma1^2/2);
    omegaIM = sqrt(pi)/(2*kappa-3)^(3/2) * gamma(kappa+1)/gamma(kappa-1/2) * 2^(3/2)/(sigma1^3*k^3)*( ...
        (1-beta)/beta * (sigma1/sigma2)^3 * (k*v0/omegaRE-1) * (1+2/((2*kappa-3)*k^2*sigma2^2) * (1-k*v0/omegaRE)^2)^(-(kappa+1))...
        -( 1+2/((2*kappa-3)*k^2*sigma1^2) + 3/(2*kappa-3) )^(-(kappa+1)) );

    omega = omegaRE + 1i*omegaIM;
end