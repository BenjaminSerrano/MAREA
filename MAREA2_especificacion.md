# MAREA-2 — Especificación de diseño y plan de ablación

*Metaheurística bacteriana por microcolonias, inspirada en el ciclo de vida de* Cobetia sp. *MM1IDA2H-1*
*Versión 0.2 · septiembre 2026*

---

## 1. Por qué MAREA-2 (lo que dicen los datos del OE1)

El benchmark de 11 algoritmos × 6 celdas × 71 funciones CEC × 51 corridas deja tres hechos:

| Hecho | Evidencia |
|---|---|
| **El reescalado del ruido explica casi toda la mejora de MAREA.** | MAREA_REL vs MAREA: 70/0/1, 69/1/1, 63/8/0, 70/0/1 (B/E/W, 5k y 50k). La formulación biológica de MAREA-bio con el ruido ya reescalado solo ayuda con 5k FES (45/11/15, 49/2/20) y perjudica con 50k–500k (5/7/59, 4/1/66). |
| **MAREA_REL converge pero no escapa.** | Al pasar de 50k a 500k FES el error en unimodales baja a 0.02–0.05×, en multimodales solo a 0.89–0.92×. ABC baja a 0.46–0.71×. 10× más presupuesto casi no sirve: la población está atrapada en un óptimo local sin mecanismo de salida. |
| **Contra ABC pierde exactamente en multimodales y composición.** | D50/500k: unimodal 3/1/3, multimodal 6/1/17, híbrida 7/1/11, composición 6/0/15. |

Causa estructural: la ecuación (1) de AGEO tiene **un solo atractor** (X_g y µ globales). Toda la población cae en la misma cuenca. ABC gana porque cada fuente de alimento es un atractor local independiente que se **abandona** tras `limit` iteraciones sin mejora.

Además, MAREA-bio falló por tres razones ya diagnosticadas: (i) B̃ satura en la iteración 2 → la Fase I no existe; (ii) el salto de Lévy sale desde X_g → destruye la convergencia del élite (unimodales 0.61–0.69× al pasar de 50k a 500k); (iii) la dispersión es probabilística (p_det), no disparada por el estado de la búsqueda.

## 2. Principios de diseño

1. **Mecanismo primero, metáfora después.** Cada componente responde a una carencia medida y luego recibe su nombre biológico. Es la única defensa sólida ante la crítica de Sörensen (2015, *Metaheuristics — the metaphor exposed*).
2. **Un flag por componente.** Nada entra al algoritmo final sin sobrevivir la ablación (sección 5).
3. **Nunca unidades absolutas.** Todo paso es relativo a la geometría de la colonia o al dominio. El bug de γ no puede repetirse.
4. **El élite nunca se destruye.** Ninguna perturbación sale desde X_g; la dispersión afecta a las células peores de una colonia estancada.
5. **Los disparadores dependen del estado de la búsqueda, no de t/T.** Estancamiento, tasa de éxito, densidad — no umbrales que se cruzan al segundo paso.

## 3. Componentes

| # | Carencia que ataca | Mecanismo de búsqueda | Análogo biológico (Cobetia) | Flag |
|---|---|---|---|---|
| C1 | Un solo atractor → colapso a una cuenca | **Microcolonias**: k-means dinámico; X_c*, µ_c, σ_c por colonia. La ec. (1) se aplica con los valores locales. | Biofilm formado por microcolonias, no una biopelícula homogénea | `colonies` |
| C2 | Sin mecanismo de escape | **Inanición → dispersión**: contador τ_c de iteraciones sin mejorar X_c*. Si τ_c > limit, la fracción p peor de la colonia salta (Lévy, desde su propia posición); la mejor célula se conserva. | Desprendimiento del biofilm disparado por agotamiento de nutrientes; células persistentes | `starvation` |
| C3 | Paso mal calibrado / colapso auto-reforzado de σ | **Ley de trade-off µ_max–biofilm**: regla 1/5 por colonia. γ_c = γ · exp(k·(s_c − 0.2)), s_c = EMA de la tasa de éxito. Éxito alto = régimen planctónico (paso mayor); éxito bajo = régimen sésil (paso menor, y si persiste, desprendimiento C2). | Trade-off crecimiento–biofilm del paper (ec. 7, Fig. 1, sec. 3.5): la célula invierte en proliferar cuando el ambiente rinde y en estructura cuando no | `quorum` |
| C4 | Población fija ineficiente (LPSR de L-SHADE es de lo poco que se sabe que funciona) | **Logística de Cobetia** como capacidad de carga: N(u) = N₀ − (N₀ − N_min)·G(u), con G la logística en log-tiempo **ajustada en el paper** (ec. 1, Tabla 1: x_mid = 3.429, s = 0.784, ventana 63 h), normalizada a G(1)=1. | Consumo de sustrato en cultivo batch. Es la única pieza que toma literalmente la forma funcional y los números del modelo GAMLSS. | `pop_schedule='cobetia'` |
| C5 | La Fase I no existía | **Fase de latencia**: durante los primeros λ·MaxFES, α=α_min, β=β_max, N=N₀. | Fase lag del crecimiento bacteriano | `lag` |
| C6 | Re-exploración de cuencas ya agotadas | **Memoria de nutrientes**: archivo de centros de colonias dispersadas; células que caen a < r de uno se re-siembran. | Zonas de sustrato agotado quedan inhabitables | `nutrient_memory` |
| C7 | Ruido por eje no sigue valles rotados (CEC está lleno de rotaciones) | **Perturbación por diferencia de vecinas**: ξ = (X_a − X_b) de dos células de la misma colonia + 0.3·σ_c·t_ν, con **t_ν Student-t de ν = 3.05**, la familia TF que el paper ajustó a los residuos de crecimiento (Tabla 1). Invariante a rotación sin estimar covarianza; colas pesadas. | Transferencia horizontal de genes entre vecinas; variabilidad de colas pesadas medida en Cobetia | `noise='diff'`, `axis_dist='student'` |
| C8 | La ec. (1) mueve sin seleccionar | **Aceptación greedy** por célula. | Selección natural | `greedy` |

### Ecuación de movimiento (célula i de la colonia c)

```
X_i ← X_i + α_t (X_c* − X_i) + β_t (µ_c − X_i) + γ_c · ξ_i
ξ_i = (X_a − X_b) · N(1, 0.3)  +  0.3 · σ_c ⊙ t_ν(I),  ν = 3.05     a ≠ b ∈ colonia c
α_t = α_min + (α_max − α_min)·v ,  β_t = β_max (1 − v) ,  γ_c = γ_max (1 − v + 0.05) · q_c
v   = 0 si u < λ ;  (u − λ)/(1 − λ) si no          (u = FES/MaxFES)
q_c = clip( exp( k_q (s_c − 0.2) ), 0.25, 4 )
```

Con `colonies=False, starvation=False, quorum=False, pop_schedule='constant', noise='uniform', greedy=False, lag=0, α_max=0.99, β_max=0.3` el algoritmo **es exactamente MAREA_REL**. Esa es la configuración `AGEO_REL` del harness y sirve de control.

## 4. Pseudocódigo

```
entrada: f, [L,U], N0, MaxFES, semilla
X ~ U(L,U)^{N0×D};  F = f(X);  fes = N0;  élite = argmin F
colonias ← kmeans(X, k = ⌈N0/m⌉);  τ_c = 0;  s_c = 0.2;  memoria = ∅
mientras fes < MaxFES:
    u = fes/MaxFES
    N_t = capacidad_de_carga(u)                              # C4 (Gompertz) / C5 (λ)
    si N_t < |X|: eliminar las |X|−N_t peores (nunca el élite); re-clusterizar
    si t mod 10 = 0: re-clusterizar (semillas = mejores de cada colonia)
    α, β, γ ← schedule(u, λ)
    para cada colonia c:
        X_c*, µ_c, σ_c ;  q_c ← regla 1/5(s_c)                # C1, C3
        X'_i ← ecuación de movimiento para i ∈ c               # C7
    X' ← reflejar en [L,U]
    F' ← f(X') (solo hasta agotar MaxFES);  fes += |X'|
    aceptar X'_i si F'_i < F_i                                # C8
    para cada colonia c:
        s_c ← 0.9 s_c + 0.1 · frac(mejoraron)
        si min F_c mejoró: τ_c = 0  si no: τ_c += 1
        si τ_c > limit:                                        # C2
            memoria ← memoria ∪ {X_c*}                         # C6
            dispersar la fracción p peor: X_i += Lévy(1.5) · 0.1(U−L)
            si X_i cae a < r de un punto de la memoria: X_i ~ U(L,U)
            evaluar dispersadas; τ_c = 0; s_c = 0.2
    actualizar élite
devolver élite
```

Presupuesto: controlado en **evaluaciones**, no en iteraciones. Con reducción de población el número de iteraciones supera T = MaxFES/N₀. El smoke test debe verificar `r['fes'] <= MaxFES`, no `r['iters'] == T`.

## 5. Parámetros por defecto (v0.1)

Fijados en un smoke set sintético (esfera, elipsoide rotado 1e6, Rastrigin, Schwefel, Griewank-Rosenbrock; D10 y D30; 50k FES). **No son definitivos**: se validan en la ablación CEC y quedan fijos para todo el protocolo.

| Parámetro | Valor | Nota |
|---|---|---|
| N₀ / N_min | 80 / 20 | N₀ lo fija el runner |
| m (células por colonia) | 20 | ⇒ 4 colonias iniciales en N₀=80 |
| λ (latencia) | 0.05 | fase lag garantizada |
| schedule: x_mid / s / t_max | 3.429 / 0.784 / 63 h | Tabla 1 del paper; G(0.5) = 0.71, G(1) = 1 |
| ν (Student-t) | 3.05 | Tabla 1 del paper |
| α_min / α_max | 0.10 / **0.50** | 0.99 (AGEO) colapsa la colonia sobre X_c* |
| β_max | **0.10** | 0.3 (AGEO) fue el mayor freno en el smoke set |
| γ_max | 0.80 | multiplica el vector diferencia |
| mezcla gaussiana en `diff` | 0.3 | cubre el subespacio no abarcado cuando |c| ≤ D |
| regla 1/5: objetivo / ganancia / clip | 0.2 / 2.0 / [0.25, 4] | |
| limit (inanición) | 10 + 2D iteraciones | D10→30, D50→110 |
| p (fracción dispersada) / escala Lévy | 0.5 / 0.1·(U−L) | |
| radio de agotamiento / tamaño memoria | 0.05 / 50 | distancia normalizada por dimensión |

Resultado del smoke set (D10, 50k FES, media de 3 semillas, error):

| config | esfera | elipsoide | Rastrigin | Schwefel | Griew-Rosen |
|---|---|---|---|---|---|
| **FULL** | 1.1e-14 | 1.7e4 | 6.0 | 587 | 0.43 |
| AGEO_REL (≈MAREA_REL) | 1.5e-1 | 5.0e4 | 31.2 | 2282 | 0.89 |

Limitación conocida: el elipsoide mal condicionado sigue lejos de lo que logra DE/CMA-ES. Sin adaptación de covarianza es esperable; queda como componente candidato futuro (C9, "covarianza por colonia"), no para v0.1.

## 6. Plan de ablación

**Fase A — screening (barato).** D10 y D50, MaxFES = 5k y 50k, 71 funciones, 15 corridas. 14 configuraciones: FULL, `sin_colonias`, `sin_inanicion`, `sin_quorum`, `sin_schedule`, `lineal` (LPSR), `gompertz`, `sin_paso_relativo`, `sin_memoria`, `sin_greedy`, `sin_lag`, `ruido_normal`, `ruido_gauss`, `AGEO_REL`. Tres de ellas contrastan directamente la fidelidad al paper: `gompertz` y `lineal` contra la logística de Cobetia, `ruido_gauss` contra la Student-t. Coste ≈ 14 × 2 × 2 × 71 × 15 ≈ 60 000 corridas, la mayoría de 5k–50k FES: una tarde de máquina.

```
python ablacion_marea2.py run    --dims 10 50 --maxfes 5000  --runs 15 --workers 8
python ablacion_marea2.py run    --dims 10 50 --maxfes 50000 --runs 15 --workers 8
python ablacion_marea2.py report --dims 10 50 --maxfes 50000
python ablacion_marea2.py ratio  --dims 10 50 --fes_a 5000 --fes_b 50000
```

**Criterio de aceptación por componente.** C_k se queda si FULL vs `sin_C_k` cumple **(a)** B > W en la clase que C_k ataca (C1, C2, C6 → multimodal + composición; C3, C7, C8 → todas; C4, C5 → híbridas y composición) **y (b)** W ≤ B + 2 en unimodales. Si C_k no cumple (a) en ninguna de las 4 celdas, se elimina. Si cumple (a) pero no (b), se retoca su parámetro y se repite solo esa comparación.

**Métrica de escape.** `ratio` = mediana err(50k)/err(5k) por clase. Un componente que ataca lo multimodal debe bajar el ratio en multimodal/composición respecto a `sin_C_k`. Con 500k disponible, repetir con `--fes_a 50000 --fes_b 500000`: el objetivo es dejar el 0.89–0.92 de MAREA_REL por debajo del 0.46–0.71 de ABC.

**Fase B — confirmación.** La configuración superviviente, D10 y D50, 500k FES, 15 corridas, contra MAREA_REL, ABC, DE y (nuevo) L-SHADE. Si no queda en el grupo de cabeza en D50, volver a Fase A antes de gastar el protocolo completo.

**Fase C — protocolo completo (setup experimental).** 51 corridas, 2 dimensiones × 4 MaxFES (incluida la celda de 5M), pool con FA/CLONALG corregidos + L-SHADE, jSO y al menos un ganador CEC reciente (EA4eig o AGSK). Friedman + Shaffer, Wilcoxon B/E/W, rankings separados por MaxFES.

## 7. Posicionamiento y novedad

Ya existen metaheurísticas bacterianas: BFO (Passino 2002, quimiotaxis + reproducción + eliminación-dispersión), Bacterial Colony Optimization (Niu & Wang 2012), Bacterial Chemotaxis (Müller 2002), y variantes. La novedad de MAREA-2 **no** es la bacteria; es el conjunto de mecanismos y su justificación empírica:

- **Dispersión disparada por inanición local** (por colonia, por estancamiento) en lugar de la eliminación-dispersión probabilística de BFO.
- **Capacidad de carga con la logística ajustada a Cobetia** (forma funcional y parámetros de la Tabla 1 del paper). Hipótesis falsable: logística de Cobetia ≥ reducción lineal (LPSR) y ≥ Gompertz genérica. Las configuraciones `lineal` y `gompertz` la contrastan directamente.
- **Perturbación Student-t con ν = 3.05**, la distribución que el paper ajustó a la variabilidad del crecimiento. Hipótesis falsable: colas pesadas ≥ gaussiana (`ruido_gauss`).
- **Regla de paso = ley de trade-off µ_max–biofilm** del paper (ec. 7), leída como ley de control: rendimiento alto → estrategia planctónica, rendimiento bajo → estrategia sésil y, si persiste, desprendimiento.
- **El estudio 2×2 del OE1** como evidencia de que la contribución biológica de MAREA-bio era un artefacto de escala. Esto es publicable por sí mismo y refuerza la credibilidad del diseño de MAREA-2: cada componente llega con su ablación.

Nicho honesto si el rendimiento en 500k no alcanza al grupo L-SHADE/ABC: **optimización costosa (≤ 5k–50k FES)**, donde la formulación biológica ya gana hoy y MAREA_REL es 2º en D50/50k.

## 8. Qué se hereda del paper y qué no (revisión de la cadena GAMLSS → AGEO → MAREA-2)

GAMLSS es un marco de **regresión** (localización, escala y forma como funciones de covariables), no un modelo dinámico. No existe una traducción directa "GAMLSS → ecuación de actualización". Lo que sí se transfiere son las **leyes ajustadas** con él:

| Del paper | Ecuación / dato | En MAREA-2 |
|---|---|---|
| Ley de crecimiento logística en log-tiempo | ec. 1, Tabla 1 (A, x_mid, s) | C4: capacidad de carga N(u) |
| Trade-off µ_max–biofilm | ec. 7, Fig. 1, sec. 3.5 | C3: regla de paso por tasa de éxito; C2: desprendimiento cuando el régimen sésil no rinde |
| Residuos de colas pesadas, TF con ν ≈ 3 | Tabla 1 | C7: forma de la perturbación |
| Dispersión no constante (depende de log t y de covariables) | ec. 3 | La escala del paso no es constante: depende de σ_c, de u y de s_c. Es el espíritu GAMLSS: modelar la escala, no solo la media. |

La ecuación (1) de AGEO (ALGO.pdf) **no** toma nada específico del paper: es un esqueleto de enjambre genérico con nombres biológicos. MAREA-2 la conserva como estructura de movimiento y le añade encima las tres leyes de arriba, que son las que hacen que el vínculo con Cobetia sea literal y falsable.

**Corrección de narrativa respecto a MAREA-bio.** MAREA-bio etiquetó la fase planctónica como explotación y la de biofilm como exploración. La biología del paper es la inversa: planctónico = células libres, proliferación rápida y dispersión (exploración); biofilm = adhesión a una superficie favorable, estructura local y supervivencia bajo estrés (explotación de una zona). MAREA-2 usa el sentido correcto: colonias = biofilm = nicho local; desprendimiento = exploración; fase lag = arranque con β alto.

## 9. Riesgos y qué hacer con ellos

| Riesgo | Mitigación |
|---|---|
| Las colonias no aportan en CEC (en el smoke set no se distinguen del ruido) | Es exactamente lo que la Fase A decide. Si `sin_colonias` ≈ FULL, C1 se cae y la dispersión pasa a ser global por estancamiento del élite. |
| Elipsoides/valles rotados mal condicionados | Reportar como limitación; C9 (covarianza por colonia, rank-µ como CMA-ES) en trabajo futuro. |
| Demasiados parámetros (≈15) | Fijarlos en Fase A y **no tocarlos** en B y C; publicar la tabla completa. La ablación muestra a cuáles es sensible. |
| CLONALG con error de nivel aleatorio incluso con presupuesto corregido | Revisar que la hipermutación decaiga con la afinidad normalizada antes de incluirlo en el paper. |

## 10. Archivos

- `marea2.py` — implementación de referencia, `MAREA2Params`, `config(nombre)`, `ABLATION_CONFIGS`.
- `ablacion_marea2.py` — harness: `run` / `report` / `ratio`; JSON incremental compatible con `experimento_config`.
- Integración en un notebook estilo `experimento_*.ipynb`:

```python
from marea2 import marea2, config
p = config('FULL')
def alg_fn(f, bounds, N, T, seed):
    p.N0 = N
    return marea2(f, bounds, p, maxfes=N*T, seed=seed)
# smoke test: assert r['fes'] <= N*T   (no r['iters'] == T)
```
