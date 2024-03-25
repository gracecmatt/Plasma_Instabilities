# README

## Notes
- Use spectral code as initial condition
- Run Xie code shifted and scaled

## Summary of Files
1. `BiKappa_AllSolScript.m` Solves the dispersion relation for all solutions and attempts to ignore redundant or incorrect solutions.
3. `BiKappa_Disp_Using_Xie.m` Solves the dispersion relation for the kappa bump-on-tail using Fourier series approximation for the complex integral as in the Xie/Weideman algorithm.
3. `BiKappa_dielectric.m` Provides the functional forms for the kappa dielectric functions for $\kappa=2,6$.
4. `BohmGross_BiKap.m` Solves the Bohm-Gross relation for this distribution (something isn't working).
1. `Global_Sens_Plots_from_Data_Instability.m` Plots the 1D sufficient summary plots with some variation from a given base parameter set. 
2. `Global_Sensitivity_PIC_VP_Kappa_dispersion_rate_p6.m` Code to run active subspaces.
2. `Global_Sensitivity_PIC_VP_Kappa_dispersion_rate_p7.m` Same as above but varying $\kappa$.
7. `slurm_GS_PIC_BiKap_p6.slurm` Use to run script on Mio.
4. `SSP_multiple_weight_vectors.m` Computes plots of 2D eigenvector projections from active subspces.
5. `TestScript_BiKap.m` Runs all of the tests needed to trust the methods and run active subspaces
5. `Vlasov_1D_linearized_Steve_v4_Kappa.m` Fourth edition of the spectral method for solving VP with kappa equilibrium distribution.
8. `zetaph.m` Code to produce Fourier series approximation as in Xie/Weideman algorithm.
