# README

## Summary of Files
1. `Global_Sens_Plots_from_Data_Instability.m` Plots the 1D sufficient summary plots with some variation from a given base parameter set. 
2. `Global_Sensitivity_PIC_VP_Kappa_dispersion_rate_p6.m` Code to run active subspaces.
3. `Kappa_Bump_Disp_Using_Xie.m` Solves the dispersion relation for the kappa bump-on-tail using Fourier series approximation for the complex integral as in the Xie/Weideman algorithm.
4. `SSP_multiple_weight_vectors.m` Computes plots of 2D eigenvector projections from active subspces.
5. `Vlasov_1D_linearized_Steve_v4_Kappa.m` Fourth edition of the spectral method for solving VP with kappa equilibrium distribution.
6. `dielectric_kappa.m` Provides the functional forms for the kappa dielectric functions for $\kappa=1,2,6$.
7. `slurm_GS_PIC_BiKap_p6.slurm` Use to run script on Mio.
8. `zetaph.m` Code to produce Fourier series approximation as in Xie/Weideman algorithm.

## Sufficient Summary Plots
### Kappa = 1
#### Bi-Kappa $(\kappa=1)$ with 1% Variation
![1% Variation](Figs/kappa1/EigWVSSPfit_Dispersion_KappaBump_1_512_2.svg)
#### Bi-Kappa $(\kappa=1)$ with 5% Variation
![5% Variation](Figs/kappa1/EigWVSSPfit_Dispersion_KappaBump_5_512_2.svg)
#### Bi-Kappa $(\kappa=1)$ with 10% Variation
![10% Variation](Figs/kappa1/EigWVSSPfit_Dispersion_KappaBump_10_512_2.svg)
#### Bi-Kappa $(\kappa=1)$ with 15% Variation
![15% Variation](Figs/kappa1/EigWVSSPfit_Dispersion_KappaBump_15_512_2.svg)
#### Bi-Kappa $(\kappa=1)$ with 25% Variation
![25% Variation](Figs/kappa1/EigWVSSPfit_Dispersion_KappaBump_25_512_2.svg)
#### Bi-Kappa $(\kappa=1)$ with 50% Variation
![50% Variation](Figs/kappa1/EigWVSSPfit_Dispersion_KappaBump_50_512_2.svg)

### Kappa = 2
#### Bi-Kappa $(\kappa=2)$ with 1% Variation
![1% Variation](Figs/kappa2/EigWVSSPfit_Dispersion_KappaBump_1_512_2.svg)
#### Bi-Kappa $(\kappa=2)$ with 5% Variation
![5% Variation](Figs/kappa2/EigWVSSPfit_Dispersion_KappaBump_5_512_2.svg)
#### Bi-Kappa $(\kappa=2)$ with 10% Variation
![10% Variation](Figs/kappa2/EigWVSSPfit_Dispersion_KappaBump_10_512_2.svg)
#### Bi-Kappa $(\kappa=2)$ with 15% Variation
![15% Variation](Figs/kappa2/EigWVSSPfit_Dispersion_KappaBump_15_512_2.svg)
#### Bi-Kappa $(\kappa=2)$ with 25% Variation
![25% Variation](Figs/kappa2/EigWVSSPfit_Dispersion_KappaBump_25_512_2.svg)
#### Bi-Kappa $(\kappa=2)$ with 50% Variation
![50% Variation](Figs/kappa2/EigWVSSPfit_Dispersion_KappaBump_50_512_2.svg)

### Kappa = 6
#### Bi-Kappa $(\kappa=6)$ with 1% Variation
![1% Variation](Figs/kappa6/EigWVSSPfit_Dispersion_KappaBump_1_512_2.svg)
#### Bi-Kappa $(\kappa=6)$ with 5% Variation
![5% Variation](Figs/kappa6/EigWVSSPfit_Dispersion_KappaBump_5_512_2.svg)
#### Bi-Kappa $(\kappa=6)$ with 10% Variation
![10% Variation](Figs/kappa6/EigWVSSPfit_Dispersion_KappaBump_10_512_2.svg)
#### Bi-Kappa $(\kappa=6)$ with 15% Variation
![15% Variation](Figs/kappa6/EigWVSSPfit_Dispersion_KappaBump_15_512_2.svg)
#### Bi-Kappa $(\kappa=6)$ with 25% Variation
![25% Variation](Figs/kappa6/EigWVSSPfit_Dispersion_KappaBump_25_512_2.svg)
#### Bi-Kappa $(\kappa=6)$ with 50% Variation
![50% Variation](Figs/kappa6/EigWVSSPfit_Dispersion_KappaBump_50_512_2.svg)
