function omega = analyticAPPROX(k,sigma,mu)

    v = -500:0.1:500;
    f0 = @(z) exp(-(z-mu).^2/sigma^2)/sqrt(pi*sigma^2);
    f0p = @(z) exp(-(z-mu).^2/sigma^2)/sqrt(pi*sigma^2).*(2*(z-mu)/sigma^2);

    % REAL PART (Omega)
    C2 = 3*trapz(v.^2.*f0(v)); %if f0 is even
    Omega = sqrt(1+C2.*k^2);

    % IMAGINARY PART (gamma)
    gamma = Omega*(-pi/2)*(1/k^2)*f0p(Omega/k);

    omega = Omega+1i*gamma;

end