# Summary of Kappa Equilibrium Distribution for Stability Analysis via Linearized Dispersion Relation

**Grace Mattingly**  
**August 1, 2023**

A 3D single-species kappa equilibrium velocity distribution function is given by:

$$\begin{align}\displaystyle 
    f_\kappa(\bm{v};\bm{\mu},\sigma) &= (\pi\sigma^2)^{-3/2} A_\kappa \left[1+\frac{|\bm{v}-\bm{\mu}|^2}{\sigma^2(\kappa-\frac{3}{2})}\right]^{-(\kappa+1)} \\
    A_\kappa &= \left(\kappa-\frac{3}{2}\right)^{-3/2} \frac{\Gamma(\kappa+1)}{\Gamma(\kappa-\frac{1}{2})}
\end{align}$$

where $\bm{v}$ is the particle velocity, $\bm{v}_0$ is the mean/bulk velocity, $\sigma=\sqrt{2k_B T/m}$ is the thermal velocity, $T$ is the temperature and second velocity moment, $k_B$ is the Boltzmann constant, $m$ is the particle mass, and $\kappa\in(3/2,\infty)$ is the index of the distribution.

Integrating out two of the dimensions, leaving $v=v_1$ and assuming $\bm{\mu}=\mu \hat{v}_1$ for simplicity.

$$\begin{align}\displaystyle 
    f_\kappa(v;\mu,\sigma)&=\frac{1}{\sqrt{\pi(\kappa-\frac{3}{2})\sigma^2}}\frac{\Gamma(\kappa)}{\Gamma(\kappa-\frac{1}{2})}\left[1+\frac{|v-\mu|^2}{(\kappa-\frac{3}{2})\sigma^2}\right]^{-\kappa} 
\end{align}$$

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

## Xie/Weideman Algorithm
The following values of $\kappa$ did **not** result in a solution:

- $1.75$
- $1.5$
- $4/3$
- $\pi$
