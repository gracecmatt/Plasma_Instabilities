function omega = BohmGross_Kap(k,sigma,mu,kappa)
% approximation works for small k
% assume mu = 0;

    A=(pi*sigma^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
    f0= @(z) A*(1+z.^2/((kappa-1.5)*sigma^2)).^(-kappa);
    f0p= @(z) A*2*z.*(-kappa)./((kappa-1.5)*sigma^2).*(1+z.^2/((kappa-1.5)*sigma^2)).^(-kappa-1);

    % REAL PART (Omega)
    C1 = 0;% -2*integral(@(v) v.*f0(v),-Inf,Inf); %second moment
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

end