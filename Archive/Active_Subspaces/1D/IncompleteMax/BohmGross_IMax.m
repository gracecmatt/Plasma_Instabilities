function omega = BohmGross_IMax(k,sigma,mu,nu)
% Maxwellian; approximation works for small k

%     f0 = @(z) 1/sqrt(pi*sigma^2)*exp(-(z-mu).^2/sigma^2);
%     f0p = @(z) 2.*(z-mu)/sigma^2 .* 1/sqrt(pi*sigma^2).*exp(-(z-mu).^2/sigma^2);
    f0 = @(z) 2*heaviside(z-nu)/sqrt(pi*sigma^2).*exp(-(z-mu).^2/(sigma^2))/(1+erf(mu-nu)/sigma);
    f0p = @(z) 2*dirac(z-nu).*exp(-((z-mu)/sigma).^2)/(sqrt(pi*sigma^2)*(1+erf((mu-nu)/sigma))) ...
        - 4*(z-mu).*heaviside(z-nu).*exp(-((z-mu)/sigma).^2)/(sqrt(pi*sigma^6)*(1+erf((mu-nu)/sigma)));
    
    % REAL PART (Omega)
    C1 = -2*integral(@(v) v.*f0(v),-Inf,Inf); %second moment
    C2 = 3*integral(@(v) v.^2.*f0(v),-Inf,Inf); %third moment
    syms x 
    Omega = double(solve(x^4-x^2+C1*k*x-C2*k^2==0,x,'Real',true));
    Omega = max(Omega);
    % Omega has conjugate parts, +/-

    % IMAGINARY PART (gamma)
    gamma = -Omega.*(-pi/2)*(1/k^2).*f0p(Omega/k);

    % SOLUTION
    omega = Omega+1i*gamma;
    
    % Sort roots by greatest imaginary part
    temp = sort(1i*omega,'ComparisonMethod','real');
    omega = -1i*temp(1);
    g=2;
end