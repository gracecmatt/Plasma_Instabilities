# README

## Notes:
- Use Bohm-Gross relation for initial guess
- Use shifted and scaled version for Xie code
- If using spectral code: change Vmax to 50 for $\sigma=1$ and 150 or 200 for $\sigma=5$.

## Summary of Files
1. `BohmGross_Kap.m` computes an approximate solution for $\omega$, one that is only accurate for small $k$.
1. `Global_Sens_Plots_from_Data_Instability.m` Plots the 1D sufficient summary plots with some variation from a given base parameter set. 
2. `Global_Sensitivity_PIC_VP_Kappa_dispersion_rate_p4.m` Code to run active subspaces with 4 parameters: $k,\theta,\mu,\kappa$.
3. `Kappa_AllSolScript.m` Solves the dispersion relation for all the roots and attempts to ignore redundant or fake roots.
3. `Kappa_dielectric.m` Provides the functional forms for the kappa dielectric functions for $\kappa=2,6$.
3. `Kappa_Disp_Using_Xie.m` Solves the dispersion relation for the kappa bump-on-tail using Fourier series approximation for the complex integral as in the Xie/Weideman algorithm.
7. `slurm_GS_PIC_BiKap_p4.slurm` Use to run script on Mio.
5. `TestScript_Kap.m` Runs all of the tests needed to trust the methods and run active subspaces.
6. `TestScript_KapVar.m` Same as above, but changes $\kappa$ as a parameter.
6. `Vlasov_1D_linearized_Steve_v4_Kappa.m` Fourth edition of the spectral method for solving VP with kappa equilibrium distribution.
8. `zetaph.m` Code to produce Fourier series approximation as in Xie/Weideman algorithm.
