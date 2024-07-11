# Incomplete Bi-Maxwellian
**Distribution**
$$f_{IBM}(v)= H(v-\nu) \left( \frac{\beta \exp\left[-\frac{(v-\mu_1)^2}{\sigma_1^2}\right]}{\sqrt{\pi \sigma_1^2} \left(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1})\right)\frac{1}{2}} + \frac{(1-\beta) \exp\left[-\frac{(v-\mu_2)^2}{\sigma_2^2}\right]}{\sqrt{\pi \sigma_2^2} \left(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2})\right) \frac{1}{2}} \right)$$

**Derivative**
$$f_{IBM}'(v)=\frac{\beta \exp\left[-\frac{(v-\mu_1)^2}{\sigma_1^2}\right] \left(\sigma_1^2\delta(v-\nu)-2(v-\mu_1)H(v-\nu)\right) }{\sqrt{\pi}\sigma_1^3 \left(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1})\right)\frac{1}{2}} \\ + \frac{(1-\beta) \exp\left[-\frac{(v-\mu_2)^2}{\sigma_2^2}\right] \left(\sigma_2^2\delta(v-\nu)-2(v-\mu_2)H(v-\nu)\right)}{\sqrt{\pi}\sigma_2^3 \left(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2})\right)\frac{1}{2}}$$

**Integral to compute using Xie/Weideman algorithm**
```math
\begin{aligned}
\int \frac{f_{IBM}'(v)}{v-z}dv&= \left. \left[ \frac{1}{v-z}\frac{\beta}{\sqrt{\pi}\sigma_1}\frac{\exp\left[\frac{-(v-\mu_1)^2}{\sigma_1^2}\right]}{\left(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1})\right)\frac{1}{2}} +  \frac{1}{v-z}\frac{(1-\beta)}{\sqrt{\pi}\sigma_2}\frac{\exp\left[\frac{-(v-\mu_2)^2}{\sigma_2^2}\right]}{\left(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2})\right)\frac{1}{2}} \right]  \right|_{v=\nu} \\ 
&+\int \frac{H(v-\nu)}{v-z} \left\{ \frac{-2\beta(v-\mu_1)}{\sqrt{\pi}\sigma_1^3}\frac{\exp\left[\frac{-(v-\mu_1)^2}{\sigma_1^2}\right]}{\left(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1})\right)\frac{1}{2}} + \frac{-2(1-\beta)(v-\mu_2)}{\sqrt{\pi}\sigma_2^3}\frac{\exp\left[\frac{-(v-\mu_2)^2}{\sigma_2^2}\right]}{\left(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2})\right)\frac{1}{2}} \right\}dv
\end{aligned}
```
Simplified
```math
\begin{aligned}
\int \frac{f_{IBM}'(v)}{v-z}dv&=\text{Zp} + \frac{2}{\nu-z} \left\{ \frac{\beta}{\sqrt{\pi}\sigma_1}\frac{\exp\left[\frac{-(\nu-\mu_1)^2}{\sigma_1^2}\right]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} +  \frac{(1-\beta)}{\sqrt{\pi}\sigma_2}\frac{\exp\left[\frac{-(\nu-\mu_2)^2}{\sigma_2^2}\right]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))} \right\}, \\ 
\text{Zp} &= \text{calZ}\left(z, 
-4 H(v-\nu) \left\{ \frac{\beta(v-\mu_1)}{\sqrt{\pi}\sigma_1^3}\frac{\exp\left[\frac{-(v-\mu_1)^2}{\sigma_1^2}\right]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} + \frac{(1-\beta)(v-\mu_2)}{\sqrt{\pi}\sigma_2^3}\frac{\exp\left[\frac{-(v-\mu_2)^2}{\sigma_2^2}\right]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\right\}, N\right)
\end{aligned}
```
