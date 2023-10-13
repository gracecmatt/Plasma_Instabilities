# README

## To Do: Friday, 10/13/23
- Update/add equations for each subfolder (kappa is wrong)
- add rescaled level curves to kappa, kappa var, incomplete maxwellian
- save level curves in Figs folder


## Sub-Folders:
- 1D BiMaxwellian $f_{BM}(v) = \frac{\beta}{\sqrt{\pi\sigma_1^2}}\text{exp}\left[-\frac{|v-\mu_1|^2}{\sigma_1^2}\right] + \frac{1-\beta}{\sqrt{\pi\sigma_2^2}}\text{exp}\left[-\frac{|v-\mu_2|^2}{\sigma_2^2}\right]$
    - See sufficient summary plots with $k=0.5, \sigma_1=\sigma_2=0.5, \mu_1=0, \mu_2=4, \beta=0.8$ [HERE.](https://github.com/gracecmatt/Plasma_Instabilities/tree/main/Active_Subspaces/1D/BiMaxwellian/README.md)
- 1D Kappa $f_{\kappa}(v) = \frac{1}{\sqrt{\pi\theta^2\left(\kappa-\frac{1}{2}\right)}}\frac{\Gamma(\kappa+1)}{\Gamma(\kappa+1/2)}\left[1+\frac{1}{\kappa-\frac{1}{2}}\frac{({v}-{\mu})^2}{\theta^2}\right]^{-\kappa-1}$
    - See sufficient summary plots with $k=0.5, \theta=1, \mu=0$ [HERE.](https://github.com/gracecmatt/Plasma_Instabilities/tree/main/Active_Subspaces/1D/Kappa/README.md)
- 1D BiKappa $f_{B\kappa}(v) = \frac{\beta}{\sqrt{\pi\theta_1^2\left(\kappa-\frac{1}{2}\right)}} \frac{\Gamma(\kappa+1)}{\Gamma(\kappa+1/2)}\left[1+\frac{1}{\kappa-\frac{1}{2}}\frac{({v}-{\mu_1})^2}{\theta_1^2}\right]^{-\kappa-1} + \frac{1-\beta}{\sqrt{\pi\theta_2^2\left(\kappa-\frac{1}{2}\right)}}   \frac{\Gamma(\kappa+1)}{\Gamma(\kappa+1/2)}\left[1+\frac{1}{\kappa-\frac{1}{2}}\frac{({v}-{\mu_2})^2}{\theta_2^2}\right]^{-\kappa-1}$
    - See sufficient summary plots with $k=0.5, \theta_1=\theta_2=1, \mu_1=0, \mu_2=4, \beta=0.9$ [HERE.](https://github.com/gracecmatt/Plasma_Instabilities/blob/main/Active_Subspaces/1D/BiKappa/README.md)
- 1D/SummersThorne $f_{\kappa ST}(v)=\frac{1}{\sqrt{\pi\theta^2\kappa^3}}\frac{\Gamma(\kappa+1)}{\Gamma(\kappa-1/2)}\left[1+\frac{v^2}{\kappa\theta^2}\right]^{-\kappa}$

