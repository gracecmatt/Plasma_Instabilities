# Summary of Kappa Equilibrium Distribution for Stability Analysis via Linearized Dispersion Relation

**Grace Mattingly**  
**August 1, 2023**

A 3D single-species kappa equilibrium velocity distribution function is given by:

$$\begin{align}\displaystyle 
    f_\kappa(\mathbf{v};\mathbf{\mu},\sigma) &= (\pi\sigma^2)^{-3/2} A_\kappa \left[1+\frac{|\mathbf{v}-\mathbf{\mu}|^2}{\sigma^2(\kappa-\frac{3}{2})}\right]^{-(\kappa+1)} \\
    A_\kappa &= \left(\kappa-\frac{3}{2}\right)^{-3/2} \frac{\Gamma(\kappa+1)}{\Gamma(\kappa-\frac{1}{2})}
\end{align}$$

where $\mathbf{v}$ is the particle velocity, $\mathbf{v}_0$ is the mean/bulk velocity, $\sigma=\sqrt{2k_B T/m}$ is the thermal velocity, $T$ is the temperature and second velocity moment, $k_B$ is the Boltzmann constant, $m$ is the particle mass, and $\kappa\in(3/2,\infty)$ is the index of the distribution.

Integrating out two of the dimensions, leaving $v=v_1$ and assuming $\mathbf{\mu}=\mu \hat{v}_1$ for simplicity.

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

We require $\kappa>3/2$ and if $\kappa\in\mathbb{N}=\{1,2,3,\dots \}$, then this integral can be solved with the Residue Theorem. Focusing on just the integral, we have:

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

This "one-solve-all approach" seeks to find a numerical approximation for integrals of the form

$$\begin{align*}
    \int_{-\infty}^\infty \frac{F(z)}{z-\xi}dz, && \text{Im}(\xi)>0 
\end{align*}$$

Where $F(z)$ is an entire function for $z\in ?$. The integral is extended to the lower half plane by deforming the contour of integration around the pole and picking up a half or full circle as:

$$\begin{align*}
    &PV \int_{-\infty}^\infty \frac{F(z)}{z-\xi}dz + i F(z), && \text{Im}(\xi)=0 \\
    &\int_{-\infty}^\infty \frac{F(z)}{z-\xi}dz +2i F(z), && \text{Im}(\xi)<0
\end{align*}$$

Since $F(z)$ is [analytic?], we can assume an expansion 

$$\begin{align*}
   [W(z)]^{-1}F(z) = \sum_{n=-\infty}^\infty a_n \rho_n(z), && z\in\mathbb{R}
\end{align*}$$

for an orthogonal basis set $\{\rho_n(z)\}$ with corresponding weight function $W(z)$ which satisfies $\int_{-\infty}^\infty W(z)\rho_n(z)\rho_m^*(z) dz = A \delta_{n,m}$.

The coefficients are given by 

$$\begin{align*}
    a_n = \frac{1}{A} \int_{-\infty}^\infty F(z)\rho_n^*(z) dz
\end{align*}$$

and the integrand above can be expressed as 

$$\begin{align*}
   \frac{F(z)}{z-\xi} = \sum_{n=-\infty}^\infty a_n \left[ W(z) \frac{\rho_n(z)}{z-\xi} \right]
\end{align*}$$

Choosing the following basis function set and weight function:

$$\begin{align*}
   \rho_n(z)=\frac{(L+iz)^n}{(L-iz)^n}, && W(z)=\frac{1}{L^2+z^2}
\end{align*}$$

and given the transformation $z=L \tan(\theta/2)$, the coefficients can be computed with a standard FFT as $(L+iz)/(L-iz)=e^{i \theta}$. 

Then, since the coefficients are independent of $z$, the integral can be computed inside the sum using the Residue Theorem,

$$\begin{align*} \displaystyle
   \int_{-\infty}^\infty\frac{F(z)}{z-\xi} dz &=  \sum_{n=-\infty}^\infty a_n \int_{-\infty}^\infty \left[\frac{1}{L^2+z^2} \frac{1}{z-\xi} \frac{(L+iz)^n}{(L-iz)^n}\right] dz\\
    &=  \sum_{n=-\infty}^\infty a_n \begin{cases} \displaystyle
        \frac{i\pi}{L}\frac{1}{L-iz}, && n=0\\ \displaystyle
        \frac{2i\pi}{L^2 + z^2}\frac{(L+iz)^n}{(L-iz)^n}, &&n>0\\
        0, && n<0
        \end{cases}.
\end{align*}$$

Thus, 

$$\begin{align*}
   \int_{-\infty}^\infty\frac{F(z)}{z-\xi} dz = \frac{i\pi a_0}{L(L-iz)} +\frac{2i\pi}{L^2 + z^2} \sum_{n=1}^\infty a_n \left(\frac{L+iz}{L-iz}\right)^n, && \text{Im}(\xi)>0
\end{align*}$$
