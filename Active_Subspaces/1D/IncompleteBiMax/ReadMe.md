# Incomplete Bi-Maxwellian
Distribution
$$f_{IBM}(v)=2 H(v-\nu)\{\frac{\beta \exp[-\frac{(v-\mu_1)^2}{\sigma_1^2}]}{\sqrt{\pi \sigma_1^2} (1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} + \frac{(1-\beta) \exp[-\frac{(v-\mu_2)^2}{\sigma_2^2}]}{\sqrt{\pi \sigma_2^2} (1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\}$$

Derivative
$$f_{IBM}'(v)=2\{\frac{\beta \exp[-\frac{(v-\mu_1)^2}{\sigma_1^2}] (\sigma_1^2\delta(v-\nu)-2(v-\mu_1)H(v-\nu)) }{\sqrt{\pi}\sigma_1^3 (1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} \\ + \frac{(1-\beta) \exp[-\frac{(v-\mu_2)^2}{\sigma_2^2}] (\sigma_2^2\delta(v-\nu)-2(v-\mu_2)H(v-\nu))}{\sqrt{\pi}\sigma_2^3 (1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\}$$

Integral to compute using Xie/Weideman algorithm
\begin{align}
\int \frac{f_{IBM}'(v)}{v-z}dv=\left[ \frac{2\beta}{\sqrt{\pi}\sigma_1}\frac{\exp[\frac{-(v-\mu_1)^2}{\sigma_1^2}]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))}\frac{1}{v-z} +  \frac{2(1-\beta)}{\sqrt{\pi}\sigma_2}\frac{\exp[\frac{-(v-\mu_2)^2}{\sigma_2^2}]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\frac{1}{v-z} \right]\|_{v=\nu} \\ 
-4\int \frac{H(v-\nu)}{v-z} \left\{ \frac{\beta(v-\mu_1)}{\sqrt{\pi}\sigma_1^3}\frac{\exp[\frac{-(v-\mu_1)^2}{\sigma_1^2}]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} + \frac{(1-\beta)(v-\mu_2)}{\sqrt{\pi}\sigma_2^3}\frac{\exp[\frac{-(v-\mu_2)^2}{\sigma_2^2}]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\right\}dv
\end{align}

Simplified
\begin{align}
\int \frac{f_{IBM}'(v)}{v-z}dv=\text{Zp} + \frac{2}{\nu-z}\left\{ \frac{\beta}{\sqrt{\pi}\sigma_1}\frac{\exp[\frac{-(\nu-\mu_1)^2}{\sigma_1^2}]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} +  \frac{(1-\beta)}{\sqrt{\pi}\sigma_2}\frac{\exp[\frac{-(\nu-\mu_2)^2}{\sigma_2^2}]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))} \right\}, \\ 
Zp = \text{calZ}\left(z, 
-4 H(v-\nu) \left\{ \frac{\beta(v-\mu_1)}{\sqrt{\pi}\sigma_1^3}\frac{\exp[\frac{-(v-\mu_1)^2}{\sigma_1^2}]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} + \frac{(1-\beta)(v-\mu_2)}{\sqrt{\pi}\sigma_2^3}\frac{\exp[\frac{-(v-\mu_2)^2}{\sigma_2^2}]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\right\}, N\right)
\end{align}
