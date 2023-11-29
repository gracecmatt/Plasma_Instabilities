# Incomplete Bi-Maxwellian
Distribution
$$f_{IBM}(v)=2 H(v-\nu)\Big\{\frac{\beta \exp[-\frac{(v-\mu_1)^2}{\sigma_1^2}]}{\sqrt{\pi \sigma_1^2} (1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} + \frac{(1-\beta) \exp[-\frac{(v-\mu_2)^2}{\sigma_2^2}]}{\sqrt{\pi \sigma_2^2} (1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\Big\}$$

Derivative
$$f_{IBM}'(v)=2\Big\{\frac{\beta \exp\Big[-\frac{(v-\mu_1)^2}{\sigma_1^2}\Big] \Big(\sigma_1^2\delta(v-\nu)-2(v-\mu_1)H(v-\nu)\Big) }{\sqrt{\pi}\sigma_1^3 (1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} \\ + \frac{(1-\beta) \exp\Big[-\frac{(v-\mu_2)^2}{\sigma_2^2}\Big] \Big(\sigma_2^2\delta(v-\nu)-2(v-\mu_2)H(v-\nu)\Big)}{\sqrt{\pi}\sigma_2^3 (1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\Big\}$$

Integral to compute using Xie/Weideman algorithm
$$\int \frac{f_{IBM}'(v)}{v-z}dv=\Big[ \frac{2\beta}{\sqrt{\pi}\sigma_1}\frac{\exp\Big[\frac{-(v-\mu_1)^2}{\sigma_1^2}\Big]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))}\frac{1}{v-z} +  \frac{2(1-\beta)}{\sqrt{\pi}\sigma_2}\frac{\exp\Big[\frac{-(v-\mu_2)^2}{\sigma_2^2}\Big]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\frac{1}{v-z} \Big]\Big|_{v=\nu} \\ 
-4\int \frac{H(v-\nu)}{v-z} \Big\{ \frac{\beta(v-\mu_1)}{\sqrt{\pi}\sigma_1^3}\frac{\exp\Big[\frac{-(v-\mu_1)^2}{\sigma_1^2}\Big]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} + \frac{(1-\beta)(v-\mu_2)}{\sqrt{\pi}\sigma_2^3}\frac{\exp\Big[\frac{-(v-\mu_2)^2}{\sigma_2^2}\Big]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\Big\}dv$$

Simplified
$$\int \frac{f_{IBM}'(v)}{v-z}dv=\text{Zp} + \frac{2}{\nu-z}\Big\{ \frac{\beta}{\sqrt{\pi}\sigma_1}\frac{\exp\Big[\frac{-(\nu-\mu_1)^2}{\sigma_1^2}\Big]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} +  \frac{(1-\beta)}{\sqrt{\pi}\sigma_2}\frac{\exp\Big[\frac{-(\nu-\mu_2)^2}{\sigma_2^2}\Big]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))} \Big\}, \\ 
Zp = \text{calZ}\Big(z, 
-4 H(v-\nu) \Big\{ \frac{\beta(v-\mu_1)}{\sqrt{\pi}\sigma_1^3}\frac{\exp\Big[\frac{-(v-\mu_1)^2}{\sigma_1^2}\Big]}{(1+\text{erf}(\frac{\mu_1-\nu}{\sigma_1}))} + \frac{(1-\beta)(v-\mu_2)}{\sqrt{\pi}\sigma_2^3}\frac{\exp\Big[\frac{-(v-\mu_2)^2}{\sigma_2^2}\Big]}{(1+\text{erf}(\frac{\mu_2-\nu}{\sigma_2}))}\Big\}, N\Big)$$
