# Plasma Instabilities Research

## Updates:
#### 5/11/2023 GM
Started moving files to re-organize everything. The three folders to be used are `Level_Curves`, `Active_Subspaces`, and `Dielectric_Functions`. 

#### 4/13/2023 GM
Added dielectric functions for $\kappa=1,2,6$ with variable `omega` and parameters `k,beta,mu1,mu2,theta1,theta2,init_guess` to `Ben's Code/Dielectric Functions`. Functions were computing using the Mathematica notebook `Disp_Rel_Kappa_04.13.23.nb`, exported to Matlab, then formatted into usable functions calling `fsolve()` to determine the complex solution `omega`. 

#### 4/10/2023 GM
Added figures generated from `KappaBump_AllSolScript.m` plotting all unique roots in level curves for every parameter of the Kappa Bump-On-Tail with base parameter set $\kappa=1$, $\beta=0.9$, $\mu_1=0$, $\mu_2=4$, $\theta_1=1$, and $\theta_2=1$. Figures displayed in markdown file `Ben's Code/Kappa_Bump_Level_Curves/all_solutions_figures.md`.

#### 3/8/2023 GM
Added Matlab file `KappaBump_AllSolScript.m` to `Ben's Code/Level Curves k vs gamma/Kappa Bump` folder.  
This file computes all the roots for the Kappa Bump-on-Tail distribution and plots the unique roots.

#### 3/6/2023 GM
Added Ben's code for computing $\gamma$ and plotting $\gamma$ vs $k$.   
Added readme with relevant figures in `Ben's Code/Level Curves k vs gamma/Kappa Bump`.
