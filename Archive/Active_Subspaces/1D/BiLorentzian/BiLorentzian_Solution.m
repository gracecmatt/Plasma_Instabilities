function omega = BiLorentzian_Solution(k,sigma1,sigma2,mu1,mu2,beta)

    %Rescaling of dielectric function
    w0 = -sigma1*k*1i;
    w1 = k*(mu2-mu1) - sigma2*k*1i;

    wdiff = w0-w1;

    %Let y = w - w0 and compute coefficients of polynomial in y
    coeffvec = [1, 2*wdiff, wdiff^2-1, -2*beta*wdiff, -beta*wdiff^2];
    
    %Solve for roots of quartic polynomial
    rts = roots(coeffvec);

    % Transform roots of y back to roots of Omega + i*gamma
    true_rts = rts + w0 + k*mu1;
    true_rts = abs(real(rts)) + 1i*imag(rts) + w0 + k*mu1; % collapse conjugate real parts to positive one

    %Sort roots by greatest imaginary part
    temp = sort(1i*true_rts,'ComparisonMethod','real');
    true_rts = -1i*temp;

    omega = true_rts(1);
    % omega = abs(real(true_rts(1))) + 1i*imag(true_rts(1));

end

