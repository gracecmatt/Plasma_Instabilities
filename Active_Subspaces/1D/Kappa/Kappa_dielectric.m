function omega = Kappa_dielectric(k, sigma, mu, kappa, init_guess)

    if(kappa == 2)
        D = @(omega) 1+(-1).*k.^(-2).*(2.*((-1).*k.*mu+omega).^2+k.^2.*sigma.^2).^(-3).*(8.* ...
          k.^2.*((-1).*k.*mu+omega).^4+24.*k.^4.*((-1).*k.*mu+omega).^2.*sigma.^2+ ...
          (sqrt(-1)*16).*2.^(1/2).*k.^5.*(k.*mu+(-1).*omega).*sigma.^3+(-6).* ...
          k.^6.*sigma.^4);
    elseif(kappa == 6)
        D = @(omega) 1+(-2/7).*(2.*((-1).*k.*mu+omega).^2+9.*k.^2.*sigma.^2).^(-7).*(448.*(( ...
          -1).*k.*mu+omega).^12+14784.*k.^2.*((-1).*k.*mu+omega).^10.*sigma.^2+ ...
          213840.*k.^4.*((-1).*k.*mu+omega).^8.*sigma.^4+1796256.*k.^6.*((-1).*k.* ...
          mu+omega).^6.*sigma.^6+10103940.*k.^8.*((-1).*k.*mu+omega).^4.*sigma.^8+ ...
          54561276.*k.^10.*((-1).*k.*mu+omega).^2.*sigma.^10+(sqrt(-1)*60466176).* ...
          2.^(1/2).*k.^11.*(k.*mu+(-1).*omega).*sigma.^11+(-40920957).*k.^12.* ...
          sigma.^12);
    else
        disp('Error: kappa is not equal to 2 or 6');
        pause;
    end

    options = optimoptions('fsolve','Display','off');
    omega = fsolve(D, init_guess, options);
end