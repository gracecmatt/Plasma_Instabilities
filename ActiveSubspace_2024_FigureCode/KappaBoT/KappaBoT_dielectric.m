function omega = KappaBoT_dielectric(k,sigma1,sigma2,mu,v0,beta,kappa,init_guess)

    if(kappa == 2)
        D = @(omega) 1+(-1).*k.^(-2).*(2.*beta.*k.^2.*(2.*omega.^2+k.^2.*sigma1.^2).^(-3).*( ...
          4.*omega.^4+12.*k.^2.*omega.^2.*sigma1.^2+(sqrt(-1)*(-8)).*2.^(1/2).* ...
          k.^3.*omega.*sigma1.^3+(-3).*k.^4.*sigma1.^4)+2.*((-1)+beta).*k.^2.*( ...
          k.^2.*sigma2.^2+2.*(omega+(-1).*k.*v0).^2).^(-3).*(3.*k.^4.*sigma2.^4+( ...
          sqrt(-1)*8).*2.^(1/2).*k.^3.*sigma2.^3.*(omega+(-1).*k.*v0)+(-12).* ...
          k.^2.*sigma2.^2.*(omega+(-1).*k.*v0).^2+(-4).*(omega+(-1).*k.*v0).^4));
    elseif(kappa == 4)
        D = @(omega) 1+(-1).*k.^(-2).*(2.*beta.*k.^2.*(2.*omega.^2+5.*k.^2.*sigma1.^2).^(-5) ...
          .*(16.*omega.^8+224.*k.^2.*omega.^6.*sigma1.^2+1400.*k.^4.*omega.^4.* ...
          sigma1.^4+7000.*k.^6.*omega.^2.*sigma1.^6+(sqrt(-1)*(-3200)).*10.^(1/2) ...
          .*k.^7.*omega.*sigma1.^7+(-4375).*k.^8.*sigma1.^8)+2.*((-1)+beta).* ...
          k.^2.*(5.*k.^2.*sigma2.^2+2.*(omega+(-1).*k.*v0).^2).^(-5).*(4375.* ...
          k.^8.*sigma2.^8+(-7000).*k.^6.*sigma2.^6.*(omega+(-1).*k.*v0).^2+(-1400) ...
          .*k.^4.*sigma2.^4.*(omega+(-1).*k.*v0).^4+(-224).*k.^2.*sigma2.^2.*( ...
          omega+(-1).*k.*v0).^6+(-16).*(omega+(-1).*k.*v0).^8+(sqrt(-1)*(-3200)).* ...
          10.^(1/2).*k.^7.*sigma2.^7.*((-1).*omega+k.*v0)));
    else
        disp('Error: kappa is not equal to 2 or 4');
        pause;
    end

    options = optimoptions('fsolve','Display','off');
    omega = fsolve(D, init_guess, options);
end