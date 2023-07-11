# Summary of Kappa Equilibrium Distribution for Stability Analysis via Linearized Dispersion Relation

**Grace Mattingly**  
**July 11, 2023**

The 1D single-species kappa equilibrium velocity distribution function is given by:

$$\begin{align}\displaystyle 
f_\kappa(v;\mu,\theta)=\frac{1}{\sqrt{\pi\theta^2(\kappa-\frac{1}{2})}}\frac{\Gamma(\kappa+1)}{\Gamma(\kappa+\frac{1}{2})}\left[1+\frac{|v-\mu|^2}{\theta^2(\kappa-\frac{1}{2})}\right]^{-(\kappa+1)}
\end{align}$$

This distribution was chosen following the convention of this paper:
[Livadiotis and McComas (2013)](https://github.com/gracecmatt/Plasma_Instabilities/blob/main/Notes/Summaries_of_Papers.md#understanding-kappa-distributions-a-toolbox-for-space-science-and-astrophysics-livadiotis-and-mccomas-2013)

It was chosen following the advice of this editorial:
[Fichtner and Lazar (2021)](https://github.com/gracecmatt/Plasma_Instabilities/blob/main/Notes/Summaries_of_Papers.md#kappa-distributions-from-observational-evidences-via-controvesial-predictions-to-a-consistent-theory-of-nonequilibrium-plasmas-fichtner-and-lazar-2021)

> Talk about which $\kappa$ values are important and why - from experimental data

The linearized dispersion relation to be solved (derived via Fourier and Laplace transforms of the Vlasov-Poisson system) is given in dimensionless form by:

$$\begin{align}\displaystyle 
D(\omega,k)=1-\frac{1}{k^2}\int_{-\infty}^{\infty} \frac{f_{eq}'(v)}{v-\omega/k} dv=0
\end{align}$$

where $\omega=\Omega+i\gamma$ and $\Omega,\gamma\in\mathbb{R}$. The integral over the real line, $F_{eq}(\omega,k)=\int_{-\infty}^{\infty} \frac{{f_{eq}}'(v)}{v-\omega/k} dv$, must be solved by analytically continuing the integrand into the complex plane. This is done following the method of Landau. Long story short, the result will be equivalent to assuming $\text{Im}(\omega)=\gamma>0$.

Find out more in this paper by [Lev Landau (1946)](https://github.com/gracecmatt/Plasma_Instabilities/blob/main/Notes/Summaries_of_Papers.md#on-the-vibrations-of-the-electrostatic-plasma-lev-landau-1946). To equate the equations, let $\omega=is$.

The dispersion relation depends heavily on the derivative of the equilibrium distribution function of interest. Differentiating $f_\kappa(v)$ above gives:

$$\begin{align}\displaystyle 
{f_\kappa}'(v;\mu,\theta)=\frac{-2(v-\mu)(\kappa+1)}{\sqrt{\pi}\theta^3(\kappa-\frac{1}{2})^{3/2}}\frac{\Gamma(\kappa+1)}{\Gamma(\kappa+\frac{1}{2})}\left[1+\frac{|v-\mu|^2}{\theta^2(\kappa-\frac{1}{2})}\right]^{-(\kappa+2)}
\end{align}$$

We require $\kappa>1/2$ and if $\kappa\in\mathbb{N}=\{1,2,3,\dots \}$, then this integral can be solved with the Residue Theorem. Focusing on just the integral, we have:

$$\begin{align}\displaystyle 
F_\kappa(\omega,k;\mu,\theta)&=\int_{-\infty}^{\infty}\frac{-2(v-\mu)(\kappa+1)}{\sqrt{\pi}\theta^3(\kappa-\frac{1}{2})^{3/2}}\frac{\Gamma(\kappa+1)}{\Gamma(\kappa+\frac{1}{2})}\left[1+\frac{|v-\mu|^2}{\theta^2(\kappa-\frac{1}{2})}\right]^{-(\kappa+2)}\frac{dv}{v-\omega/k}\\
&\propto \int_{-\infty}^{\infty} \frac{v-\mu}{v-\omega/k}\left[\frac{1}{|v-\mu|^2+\theta^2(\kappa-\frac{1}{2})}\right]^{\kappa+2}dv
\end{align}$$

where the proportionality constant is $A(\kappa,\theta)=\frac{-2(\kappa+1)}{\sqrt{\pi}\theta^3(\kappa-\frac{1}{2})^{3/2}}\frac{\Gamma(\kappa+1)}{\Gamma(\kappa+\frac{1}{2})}\left[\theta^2\left(\kappa-\frac{1}{2}\right)\right]^{\kappa+2}$.

## Integer Values, $\kappa\in\mathbb{N}$
Using the Residue Theorem requires knowing the poles of the integrand. These are:

$$\begin{align}
p_1&=\omega/k && \text{principle root}\\
p_{2,3}&=\mu\pm i\theta\sqrt{\kappa-1/2} && \text{multiplicity $(\kappa+2)$}
\end{align}$$

We can use Mathematica's `Residue[]` function to evaluate the integrand at the poles above the Landau contour in $\mathbb{C}$. These will always be $p_1$ and $p_2$ regardless of the sign of $\text{Im}(\omega)$.

The `Integrate[]` function from Mathematica works to compute ${F}_\kappa(\omega,k)$ as well. This has been verified analytically in Mathematica by subtracting the two results and simplifying.

## Current Simulations
So far I can solve the integral $F_\kappa(\omega,k)$ for $\kappa$ values of:

- $[1,2,3,...]$

See [this folder](https://github.com/gracecmatt/Plasma_Instabilities/tree/2bceb031a948358c7603dd08ba5d2846229c07ee/Dielectric_Functions/Kappa1D) for the Mathematica-generated Matlab files of the corresponding dielectric functions.

The following values of $\kappa$ did **not** result in a solution:

- $1.75$
- $1.5$
- $4/3$
- $\pi$
