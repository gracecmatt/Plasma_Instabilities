# README

## To Do: Friday, 10/13/23
- Update/add equations for each subfolder (kappa is wrong)
- add rescaled level curves to kappa, kappa var, incomplete maxwellian
- save level curves in Figs folder

Inconsistencies:
- use sigma for all thermal velcoties
- initial guess for spectral code, use:
% Initial Condition for pertubation df(v)
omega_guess=1-1i;
df=f0_v./(omega_guess-v*k);

## Sub-Folders:
- 1D BiMaxwellian $f_{BM}(v) = \frac{\beta}{\sqrt{\pi\sigma_1^2}}\text{exp}\left[-\frac{(v-\mu_1)^2}{\sigma_1^2}\right] + \frac{1-\beta}{\sqrt{\pi\sigma_2^2}}\text{exp}\left[-\frac{(v-\mu_2)^2}{\sigma_2^2}\right]$
- 1D Kappa $f_{\kappa}(v) = \frac{1}{\sqrt{\pi\sigma^2\left(\kappa-\frac{3}{2}\right)}}\frac{\Gamma(\kappa)}{\Gamma(\kappa-1/2)}\left[1+\frac{1}{\kappa-\frac{3}{2}}\frac{({v}-{\mu})^2}{\sigma^2}\right]^{-\kappa}$
- 1D BiKappa $f_{B\kappa}(v) = \frac{\beta}{\sqrt{\pi\sigma_1^2\left(\kappa-\frac{3}{2}\right)}}\frac{\Gamma(\kappa)}{\Gamma(\kappa-1/2)}\left[1+\frac{1}{\kappa-\frac{3}{2}}\frac{({v}-{\mu_1})^2}{\sigma_1^2}\right]^{-\kappa} + \frac{1-\beta}{\sqrt{\pi\sigma_2^2\left(\kappa-\frac{3}{2}\right)}}\frac{\Gamma(\kappa)}{\Gamma(\kappa-1/2)}\left[1+\frac{1}{\kappa-\frac{3}{2}}\frac{({v}-{\mu_2})^2}{\sigma_2^2}\right]^{-\kappa}$

