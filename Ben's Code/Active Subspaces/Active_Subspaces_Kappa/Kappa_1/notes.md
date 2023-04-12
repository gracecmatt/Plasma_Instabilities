## 4/12/23
We have multiple ways to find the value of $\gamma$ for the Kappa distribution.
The first time we ran Kappa Active Subspaces, it was with spectral method v3 as an initial guess and Kappa_Disp_Relation.m to do root finding. The results are very very good.
Recently, I've been running it with spectral method v4 and Kappa_Disp_Using_Xie.m to do root finding. The results have been much less good.
I thought Kappa_Disp_Relation.m and Kappa_Disp_Using_Xie.m produced the same results, but apparently they do not. Why do we lose accuracy with the Xie code?
