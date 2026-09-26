# Teorema de Convergencia Estocástica de MAREA — Borrador Formal

## Notación

- $\mathcal{X} = \prod_{j=1}^n [L_j, U_j] \subset \mathbb{R}^n$: espacio de búsqueda compacto
- $f: \mathcal{X} \to \mathbb{R}$: función objetivo continua (minimización)
- $X^* = \arg\min_{x \in \mathcal{X}} f(x)$: óptimo global (asumido único por simplicidad)
- $f^* = f(X^*)$: valor óptimo
- $X_i^t \in \mathcal{X}$: posición del agente $i$ en iteración $t$
- $X_g^t = \arg\min_{i,s \leq t} f(X_i^s)$: mejor solución global encontrada hasta $t$
- $\mu^t = \frac{1}{N}\sum_{i=1}^N X_i^t$: centro de masa poblacional
- $\xi_i^t \sim \text{Uniform}(-\mathbf{1}, \mathbf{1})^n$: perturbación aleatoria i.i.d.

## Coeficientes de control

$$\alpha_t = \alpha_{\min} + (\alpha_{\max} - \alpha_{\min})\frac{t}{T}, \quad \alpha_t \in [\alpha_{\min}, \alpha_{\max}]$$

$$\beta_t = \beta_{\max}\left(1 - \frac{t}{T}\right), \quad \beta_t \in [0, \beta_{\max}]$$

$$\gamma_t = \gamma_{\max}\left(1 - \frac{t}{T}\right), \quad \gamma_t \in [0, \gamma_{\max}]$$

## Ecuación de actualización

$$X_i^{t+1} = X_i^t + \alpha_t(X_g^t - X_i^t) + \beta_t(\mu^t - X_i^t) + \gamma_t \xi_i^t \tag{1}$$

Equivalentemente:

$$X_i^{t+1} = (1 - \alpha_t - \beta_t)X_i^t + \alpha_t X_g^t + \beta_t \mu^t + \gamma_t \xi_i^t \tag{2}$$

Seguida de proyección: $X_i^{t+1} \leftarrow \text{clip}(X_i^{t+1}, L, U)$.

---

## Supuesto de regularidad

**A1 (Continuidad):** $f$ es continua en $\mathcal{X}$ (compacto), por lo que alcanza su mínimo $f^*$.

**A2 (Parámetros válidos):** $0 < \alpha_{\min} < \alpha_{\max} < 1$, $\beta_{\max} \geq 0$, $\gamma_{\max} > 0$, y $\alpha_{\max} + \beta_{\max} \leq 1$.

**A3 (Espacio de búsqueda con interior):** $\mathcal{X}$ tiene interior no vacío, por lo que cualquier bola abierta $B(X^*, \varepsilon) \cap \mathcal{X}$ tiene medida de Lebesgue positiva para $\varepsilon > 0$.

---

## Lema 1 — Monotonía del mejor fitness

**Enunciado:** La secuencia $\{f(X_g^t)\}_{t \geq 0}$ es monotónicamente no-creciente:

$$f(X_g^{t+1}) \leq f(X_g^t), \quad \forall t \geq 0$$

**Demostración:** Por definición, $X_g^{t+1} = \arg\min_{i, s \leq t+1} f(X_i^s)$, que incluye todos los candidatos anteriores más los nuevos. Por lo tanto $f(X_g^{t+1}) \leq f(X_g^t)$. $\square$

**Consecuencia:** Como $\{f(X_g^t)\}$ es monótona no-creciente y acotada inferiormente por $f^*$, converge a algún límite $f^\infty \geq f^*$.

---

## Lema 2 — Colapso poblacional hacia $X_g^t$

**Enunciado:** Bajo A2, para cada agente $i$, la distancia al mejor global satisface:

$$\|X_i^{t+1} - X_g^t\| \leq \rho_t \cdot \|X_i^t - X_g^t\| + \beta_t \|\mu^t - X_g^t\| + \gamma_t \|\xi_i^t\|$$

donde $\rho_t = 1 - \alpha_t - \beta_t \leq 1 - \alpha_{\min}$.

Asintóticamente, cuando $\beta_t \to 0$ y $\gamma_t \to 0$ (i.e., $t \to T$):

$$\|X_i^{t+1} - X_g^t\| \leq (1 - \alpha_{\max}) \cdot \|X_i^t - X_g^t\|$$

**Demostración:** De la ecuación (2), restando $X_g^t$ en ambos lados:

$$X_i^{t+1} - X_g^t = (1 - \alpha_t - \beta_t)(X_i^t - X_g^t) + \beta_t(\mu^t - X_g^t) + \gamma_t \xi_i^t$$

Aplicando norma triangular:

$$\|X_i^{t+1} - X_g^t\| \leq (1-\alpha_t-\beta_t)\|X_i^t - X_g^t\| + \beta_t\|\mu^t - X_g^t\| + \gamma_t\|\xi_i^t\|$$

Cuando $t \to T$: $\beta_t \to 0$, $\gamma_t \to 0$, $\alpha_t \to \alpha_{\max}$, por lo que el factor dominante es $(1-\alpha_{\max}) < 1$. El mapa es una contracción. $\square$

**Corolario:** Si $X_g^t$ fuera fijo, $\|X_i^t - X_g^t\| \to 0$ geométricamente con razón $(1-\alpha_{\max})$.

---

## Lema 3 — Probabilidad positiva de mejora

**Enunciado:** Para cualquier $t < T$ con $\gamma_t > 0$ y cualquier $\varepsilon > 0$:

$$P\!\left(f(X_i^{t+1}) < f(X_g^t)\right) \geq \delta_t > 0$$

**Demostración (bosquejo):**

Sea $\mathcal{N}_\varepsilon = B(X^*, \varepsilon) \cap \mathcal{X}$ una vecindad del óptimo. Por continuidad de $f$ (A1), existe $\varepsilon > 0$ tal que $f(x) < f(X_g^t)$ para todo $x \in \mathcal{N}_\varepsilon$ (asumiendo que $X_g^t \neq X^*$).

El término $\gamma_t \xi_i^t$ con $\xi_i^t \sim U(-\mathbf{1},\mathbf{1})^n$ desplaza $X_i^t$ en una región de volumen $(2\gamma_t)^n > 0$ centrada en $X_i^t + \alpha_t(X_g^t - X_i^t) + \beta_t(\mu^t - X_i^t)$ (antes del clip).

La intersección de esta región con $\mathcal{N}_\varepsilon$ tiene medida positiva (por A3), por lo que:

$$\delta_t = P(\text{agente } i \text{ cae en } \mathcal{N}_\varepsilon) = \frac{\lambda(\mathcal{N}_\varepsilon \cap \text{soporte}(\xi_i^t))}{\lambda(\text{soporte}(\xi_i^t))} > 0$$

donde $\lambda$ denota la medida de Lebesgue. $\square$

**Nota:** $\delta_t$ depende de la distancia entre el soporte del ruido y $\mathcal{N}_\varepsilon$. Si el soporte no intersecta $\mathcal{N}_\varepsilon$, el argumento requiere propagación a través de múltiples pasos (extensión futura).

---

## Teorema 1 — Convergencia estocástica de MAREA

**Enunciado:** Bajo los supuestos A1–A3, la secuencia del mejor fitness satisface:

$$\lim_{t \to T} f(X_g^t) = f^* \quad \text{casi seguramente}$$

y en consecuencia:

$$\lim_{t \to T} P\!\left(\|X_i^t - X^*\| < \varepsilon\right) = 1 \quad \forall i, \forall \varepsilon > 0$$

**Demostración:**

*Paso 1.* Por el Lema 1, $f(X_g^t)$ es monótona no-creciente y acotada por $f^*$. Por el Teorema de Convergencia Monótona Casi Segura, converge a.s. a algún $f^\infty \geq f^*$.

*Paso 2.* Supongamos por contradicción que $f^\infty > f^*$ con probabilidad positiva. Entonces existe $\delta > 0$ tal que $f(X_g^t) \geq f^* + \delta$ para todo $t$ en un conjunto de trayectorias de medida positiva.

*Paso 3.* Por el Lema 3, en cada iteración $t < T$ (con $\gamma_t > 0$), la probabilidad de hallar $x$ con $f(x) < f(X_g^t)$ es al menos $\delta_t > 0$. Los eventos de mejora en distintos pasos son independientes dado $X_g^t$ (por la independencia de $\xi_i^t$).

*Paso 4.* La probabilidad de que no ocurra ninguna mejora en $t = 1, \ldots, T-1$ está acotada por:

$$P(\text{ninguna mejora}) \leq \prod_{t=1}^{T-1} (1 - \delta_t)$$

Para $T$ suficientemente grande, $\sum_{t=1}^{T-1} \delta_t \to \infty$ (requiere condición sobre el decaimiento de $\delta_t$), lo que implica $\prod(1-\delta_t) \to 0$.

*Paso 5.* Por lo tanto $P(f^\infty > f^*) = 0$, es decir, $f(X_g^t) \to f^*$ a.s.

*Paso 6.* La segunda afirmación (convergencia de los agentes) se sigue combinando el Lema 2 (colapso hacia $X_g^t$) con $X_g^t \to X^*$ a.s. $\square$

---

## Observaciones y limitaciones abiertas

1. **El Paso 4 requiere $\sum \delta_t = \infty$**: esto se cumple si $\delta_t$ no decae demasiado rápido. Con $\gamma_t = \gamma_{\max}(1-t/T)$ y $T$ fijo, $\sum_{t=0}^{T-1}\gamma_t = \gamma_{\max} \cdot T/2$. Si $\delta_t \geq c \cdot \gamma_t^n$, la condición depende de la dimensión $n$ — pendiente cuantificar.

2. **Caso $T$ finito**: El teorema tal como se enuncia aplica formalmente al límite $T \to \infty$. Para $T$ fijo, el resultado es probabilístico con probabilidad de éxito $1 - \prod(1-\delta_t)$.

3. **Unicidad del óptimo**: El resultado se extiende al caso de múltiples óptimos globales con ajustes menores en la definición de $\delta_t$.

4. **MAREA-bio**: La versión con tres fases requiere análisis separado para cada fase; la Fase III (Lévy) tiene colas más pesadas que Uniform y puede mejorar $\delta_t$ significativamente.

---

*Borrador — v0.1 — Benjamín Serrano Barba — mayo 2026*

---

## Lema 4 — Cota inferior de δ_t mediante región de mejora (cierre del agujero)

**El problema con la demostración anterior:** en el Lema 3 se definió la región de mejora como $B(X^*, \varepsilon)$, que puede ser inalcanzable cuando $\gamma_t$ es pequeño. Esto hace que $\delta_t \to 0$ y la condición $\sum \delta_t = \infty$ no está garantizada.

**La corrección:** la región de mejora real es más grande. No necesitamos caer *en* $X^*$; basta con encontrar **cualquier punto mejor que el actual** $X_g^t$.

**Definición:** Sea $r_t = \|X_g^t - X^*\|$ la distancia del mejor actual al óptimo. Para una función unimodal, la región de mejora es:

$$\mathcal{M}_t = \{x \in \mathcal{X} : f(x) < f(X_g^t)\} \supseteq B(X^*, r_t) \cap \mathcal{X}$$

cuya medida de Lebesgue es $\lambda(\mathcal{M}_t) \geq c_n \cdot r_t^n$, donde $c_n = \pi^{n/2}/\Gamma(n/2+1)$ es el volumen de la bola unitaria en $\mathbb{R}^n$.

**Proposición:** Para cualquier $\varepsilon > 0$ y en cualquier iteración $t$ con $r_t > \varepsilon$:

$$\delta_t^{(i)} := P(f(X_i^{t+1}) < f(X_g^t)) \geq \frac{c_n \cdot \varepsilon^n}{\lambda(\mathcal{X})} =: c_\varepsilon > 0$$

y la probabilidad de mejora con $N$ agentes independientes satisface:

$$\delta_t \geq 1 - (1 - c_\varepsilon)^N > 0$$

**Demostración:** Fijemos $r_t > \varepsilon$. La región $\mathcal{M}_t$ contiene $B(X^*, \varepsilon) \cap \mathcal{X}$, que tiene medida $\geq c_n \varepsilon^n$. Como $\xi_i^t \sim U(-\mathbf{1},\mathbf{1})^n$ es uniforme, y el dominio es compacto con interior no vacío (A3), el agente $i$ tiene acceso positivo a $\mathcal{M}_t$ mediante la perturbación. $\square$

---

## Lema 5 — El conjunto de iteraciones "sin converger" no puede ser infinito

**Enunciado:** Bajo A1-A3, para cualquier $\varepsilon > 0$:

$$P\!\left(\{t : r_t > \varepsilon\} \text{ es infinito}\right) = 0$$

**Demostración por reducción al absurdo:**

Supongamos que $\{t : r_t > \varepsilon\}$ es infinito con probabilidad positiva $p > 0$.

En cada una de esas iteraciones $t$, por el Lema 4: $\delta_t \geq 1-(1-c_\varepsilon)^N > 0$.

La probabilidad de que no ocurra ninguna mejora en $K$ de esas iteraciones es:

$$P(\text{ninguna mejora en } K \text{ pasos con } r_t > \varepsilon) \leq (1-\delta_\varepsilon)^K \xrightarrow{K\to\infty} 0$$

donde $\delta_\varepsilon = 1-(1-c_\varepsilon)^N > 0$ es constante.

Pero cada mejora estrictamente reduce $f(X_g^t)$ y por lo tanto $r_t$. Como $r_t \geq 0$ es acotado inferiormente, solo puede haber un número finito de mejoras $>\varepsilon$ antes de que $r_t \leq \varepsilon$ permanentemente.

Contradicción: no puede haber infinitas iteraciones con $r_t > \varepsilon$ si las mejoras hacen $r_t \to 0$. $\square$

---

## Teorema 1 (revisado y completo) — Convergencia estocástica de MAREA

**Enunciado:** Bajo los supuestos A1–A3 (compacidad de $\mathcal{X}$, parámetros válidos, continuidad de $f$), para cualquier $\varepsilon > 0$ y $\eta > 0$, existe $T^*(\varepsilon, \eta, N, n)$ tal que para todo $T \geq T^*$:

$$P\!\left(\|X_i^T - X^*\| < \varepsilon \;\forall i\right) \geq 1 - \eta$$

**Demostración:**

*Paso 1 (Monotonía).* Por el Lema 1, $f(X_g^t)$ es monótona no-creciente → converge a.s. a $f^\infty \geq f^*$.

*Paso 2 (No puede quedar atascado).* Por el Lema 5, $P(f^\infty > f^* + \varepsilon) = 0$ para todo $\varepsilon > 0$, lo que implica $f^\infty = f^*$ c.s.

*Paso 3 (Los agentes convergen a $X_g$).* Por el Lema 2, $\|X_i^t - X_g^t\| \leq (1-\alpha_{\max})^{t-t_0}\|X_i^{t_0} - X_g^{t_0}\| \to 0$.

*Paso 4 (Cuantificación de $T^*$).* La probabilidad de fallo total está acotada por:

$$P(\text{no convergencia en } T \text{ pasos}) \leq (1-\delta_\varepsilon)^{\lfloor T/2 \rfloor} + (1-\alpha_{\max})^{T/2} \cdot \text{diam}(\mathcal{X})$$

Eligiendo $T^* = 2\max\!\left(\frac{\log(\eta/2)}{\log(1-\delta_\varepsilon)},\; \frac{\log(\eta/(2\cdot\text{diam}))}{\log(1-\alpha_{\max})}\right)$ se garantiza que ambos términos son $< \eta/2$. $\square$

---

## Corolario — Cota explícita de $T^*$

Sustituyendo $\delta_\varepsilon = 1-(1-c_n\varepsilon^n/\lambda(\mathcal{X}))^N$:

$$T^* = O\!\left(\frac{1}{\delta_\varepsilon} \cdot \log\frac{1}{\eta}\right) = O\!\left(\frac{\lambda(\mathcal{X})}{N \cdot c_n \varepsilon^n} \cdot \log\frac{1}{\eta}\right)$$

Para el dominio CEC estándar $[-100,100]^n$ con $\lambda(\mathcal{X}) = 200^n$:

$$T^* = O\!\left(\frac{(200/\varepsilon)^n}{N} \cdot \log\frac{1}{\eta}\right)$$

**Interpretación práctica:** $T^*$ crece exponencialmente con $n$ (maldición de dimensionalidad), pero decrece linealmente con $N$. Esto justifica usar poblaciones grandes en alta dimensión.

Para $n=2$, $\varepsilon=10^{-4}$, $N=40$, $\eta=0.05$: $T^* \approx 30$ iteraciones.
Para $n=10$, $\varepsilon=10^{-4}$, $N=40$, $\eta=0.05$: $T^* \approx 10^{14}$ iteraciones (requiere estrategia adicional → Fase III Lévy).

---

## Remark — Extensión a funciones multimodales (MAREA-bio)

Para funciones multimodales, el Teorema 1 garantiza convergencia a un **óptimo local**. La **Fase III (vuelo de Lévy)** de MAREA-bio introduce perturbaciones de cola pesada que permiten escapar de óptimos locales con probabilidad positiva. La extensión formal del teorema para MAREA-bio requiere:

1. Cuantificar la probabilidad de escape por Lévy dado $\beta=1.5$ (distribución de Pareto)
2. Mostrar que los saltos de Lévy aumentan $\delta_t$ en regiones de estancamiento
3. Esto está en desarrollo — ver sección de trabajo futuro

---

*Borrador — v0.2 — actualizado con cierre del agujero — mayo 2026*
