function omega = analyticAPPROX(k,sigma1,sigma2,mu1,mu2,beta)
% approximation works for small k
% note: 

    f0 = @(z) beta*exp(-(z-mu1).^2/sigma1^2)/sqrt(pi*sigma1^2)+(1-beta)*exp(-(z-mu2).^2/sigma2^2)/sqrt(pi*sigma2^2);
    f0p = @(z) beta*exp(-(z-mu1).^2/sigma1^2)/sqrt(pi*sigma1^2).*(2*(z-mu1)/sigma1^2)+...
           (1-beta)*exp(-(z-mu2).^2/sigma2^2)/sqrt(pi*sigma2^2).*(2*(z-mu2)/sigma2^2);

    % REAL PART (Omega)
    C1 = -2*integral(@(v) v.*f0(v),-Inf,Inf);
    C2 = 3*integral(@(v) v.^2.*f0(v),-Inf,Inf);
    syms x 
    Omega = double(solve(x^4-x^2+C1*k*x-C2*k^2==0,x,'Real',true));
    Omega = min(Omega);

    % IMAGINARY PART (gamma)
    gamma = Omega.*(-pi/2)*(1/k^2)*f0p(Omega/k);

    omega = Omega+1i*gamma;

end