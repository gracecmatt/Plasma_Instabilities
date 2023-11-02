function omega = BohmGross_BiKap(k,sigma1,sigma2,mu1,mu2,beta,kappa)
% approximation works for small k

    A1=(pi*sigma1^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
    A2=(pi*sigma2^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
    f0= @(z) beta*A1*(1+(z-mu1).^2/((kappa-1.5)*sigma1^2)).^(-kappa) + ...
         (1-beta)*A2*(1+(z-mu2).^2/((kappa-1.5)*sigma2^2)).^(-kappa);
    f0p= @(z) beta*A1*2*(z-mu1).*(-kappa)./((kappa-1.5)*sigma1^2).*(1+(z-mu1).^2/((kappa-1.5)*sigma1^2)).^(-kappa-1) + ...
          (1-beta)*A2*2*(z-mu2).*(-kappa)./((kappa-1.5)*sigma2^2).*(1+(z-mu2).^2/((kappa-1.5)*sigma2^2)).^(-kappa-1);

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

end