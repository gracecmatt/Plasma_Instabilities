# README

## To Do: Friday, 10/13/23
- add rescaled level curves to Kappa, Lorentzian, BiLorentzian, BiMaxwellian, Incomplete Maxwellian
- save level curves in Figs folder

Inconsistencies:
- use sigma for all thermal velcoties
- `num2str()` precision of 16
- input of $\mu$ (or $\mu_1$) into spectral code should always be zero
- add condition number calculation after eigs:
```m
% Compute the condition number
cond = evalues(1)/sum(evalues);
  ```
- make sure dispersion relation is using correct version of kappa
- initial guess for spectral code, use:
```m
% Initial Condition for pertubation df(v)
% Gauss
% sigma_IC=1;%Vmax/10;
% df=f_gauss(v/sigma_IC)/sigma_IC; % gives best accuracy compared to the other options
omega_guess=1-1i;
df=f0_v./(omega_guess-v*k);
```
- noise to let you know your code finished lol `load train, sound(y,Fs)`

## Sub-Folders:
- 1D BiMaxwellian $f_{BM}(v) = \frac{\beta}{\sqrt{\pi\sigma_1^2}}\text{exp}\left[-\frac{(v-\mu_1)^2}{\sigma_1^2}\right] + \frac{1-\beta}{\sqrt{\pi\sigma_2^2}}\text{exp}\left[-\frac{(v-\mu_2)^2}{\sigma_2^2}\right]$
- 1D Kappa $f_{\kappa}(v) = \frac{1}{\sqrt{\pi\sigma^2\left(\kappa-\frac{3}{2}\right)}}\frac{\Gamma(\kappa)}{\Gamma(\kappa-1/2)}\left[1+\frac{1}{\kappa-\frac{3}{2}}\frac{({v}-{\mu})^2}{\sigma^2}\right]^{-\kappa}$
- 1D BiKappa $f_{B\kappa}(v) = \frac{\beta}{\sqrt{\pi\sigma_1^2\left(\kappa-\frac{3}{2}\right)}}\frac{\Gamma(\kappa)}{\Gamma(\kappa-1/2)}\left[1+\frac{1}{\kappa-\frac{3}{2}}\frac{({v}-{\mu_1})^2}{\sigma_1^2}\right]^{-\kappa} + \frac{1-\beta}{\sqrt{\pi\sigma_2^2\left(\kappa-\frac{3}{2}\right)}}\frac{\Gamma(\kappa)}{\Gamma(\kappa-1/2)}\left[1+\frac{1}{\kappa-\frac{3}{2}}\frac{({v}-{\mu_2})^2}{\sigma_2^2}\right]^{-\kappa}$
- 1D Lorentzian $f_{L}(v)=\frac{\sigma}{\pi}\frac{1}{(v-\mu)^2+\sigma^2}$
- 1D BiLorentzian $f_{BL}(v)=\beta\frac{\sigma_1}{\pi}\frac{1}{(v-\mu_1)^2+\sigma_1^2}+(1-\beta)\frac{\sigma_2}{\pi}\frac{1}{(v-\mu_2)^2+\sigma_2^2}$
- 1D Incomplete Maxwellian $f_{IM}(v) = \frac{1}{\sqrt{\pi\sigma^2}} H(v-\nu) \text{exp}\left[-\frac{(v-\mu)^2}{\sigma^2}\right] $

## Rescaled Distributions:
- 1D BiMaxwellian $f_{BM}(v) = \frac{\beta}{\sqrt{\pi}}\text{exp}\left[-v^2\right] + \frac{1-\beta}{\sqrt{\pi(\sigma_2/\sigma_1)^2}}\text{exp}\left[-\frac{(v-(\mu_2-\mu_1)/\sigma_1)^2}{(\sigma_2/\sigma_1)^2}\right]$
