# Plan de implementación · MAREA-2 → paper OE1

**Proyecto:** MAREA (tesis doctoral UV, grupo CMBi) · **Fecha:** 25 de septiembre de 2026 · **Versión:** 4 (26-09-2026; corrige el presupuesto efectivo)
**Base:** memoria del proyecto 2026-09-25 + contraste contra `Setup Experimental para Evaluación.docx` + revisión del estado real del código (`competidores.py`, `experimento_config.py`, `marea2.py`).

---

## Resumen

El algoritmo está cerrado en sus dos leyes centrales: MAREA-2 sin_colonias fue elegida a priori (ablación A3), sobrevive a la ablación en logística en log-tiempo y ruido Student-t, y es 1º de 12 en 5 de 6 celdas. Lo que falta se divide en tres tipos:

1. **Sólo análisis, sin máquina** — media vs mediana, empate <1e-8, rankings por benchmark, Shaffer, SD de rangos, W/T/L completos, tablas §8, figuras con MAREA-2.
2. **Con máquina, barato** — confirmación de la tercera ley (`quorum`, regla 1/5 del paso, C3 de la especificación). Ya hay evidencia preliminar a favor en la Fase A (`sin_quorum` vs FULL, 15 corridas a 5k/50k: quitar el quorum empeora en las 4 celdas), pero está medida a 5k/50k y contra FULL. Falta confirmarla a 500k y sobre `sin_colonias`.
3. **Implementación + máquina** — competidores de estado del arte. **No existen en el repo**: `competidores.py` sólo tiene PSO, DE, CLONALG, GWO, WOA, ABC y FA. Hay que implementar y validar antes de usar el mini. Es lo que fija el calendario.

Los carriles corren en paralelo: el Mac mini genera; el MacBook implementa y analiza.

---

## Contraste setup vs. estado actual

| § Setup | Exige | Tienes | Estado |
|---|---|---|---|
| 1 Benchmark | 72 funciones + desagregado por benchmark | 71 (F2 CEC2017 retirada); sólo agregado. CEC2022 va en D20 dentro del conjunto "D50" (`experimento_config.py:46`), como define el §1 | Declarar la retirada y la mezcla D50/D20; desagregar (CEC2022 por benchmark = D20) |
| 2 MaxFES | 4 presupuestos (incl. 5M), 8 celdas | 3 presupuestos, 6 celdas | Divergencia declarable. 5M = semanas de mini |
| 3 Competidores | 20–25 diversos, con ganadores CEC recientes | 7 clásicos + 4 variantes propias; con 4 nuevos serían 11 | **Brecha seria.** Implementar L-SHADE/jSO; justificar 11 en el texto |
| 4 Parámetros | Valores originales; tabla de poblaciones; config única del propuesto | `POP_SIZE = 80` y `T = maxfes // N` **para todos**; FA con gamma reescalado; config única de MAREA-2 ✓ | Desviación del §4 en todos los clásicos: declararla. Para L-SHADE/jSO no es admisible: deben recibir MaxFES y derivar su población |
| 5 Protocolo | 51 corridas; métrica = **media** | 51 ✓; análisis en **mediana** | Cambiar a media (con SC los puestos no cambian) |
| 6 Rankings | Empate <1e-8; rankings por benchmark + agregados | Empate no implementado; sólo agregados | Ambos baratos |
| 7 Estadística | Friedman + Shaffer; Wilcoxon **signed-rank** W/T/L vs cada competidor; SD de rangos | Sin Shaffer; Wilcoxon **rank-sum** vs 7 de 11; sin SD. `seed = run` empareja sólo de nombre | Rank-sum por función para W/T/L; signed-rank pareado por problema (71 medias) como test global. Reportar ambos |
| 8 Presentación | Tablas por benchmark, resumen, barras por MaxFES, tablas Wilcoxon y Friedman, suplementario | Tabla 2.1, fig2 sin MAREA-2, W/T/L parcial | Datos existen; scripts fijos en 11 algoritmos |
| 9 Discusión | Dónde destaca/no; componentes por MaxFES; robustez de posición | Material ✓ | Calcular varianza de posición (SC = {1,1,1,3,1,1}) |
| 10 Reproducibilidad | Código, semillas, versiones de competidores, crudos | JSON ✓; opfunu 1.0.4 ✓; **no hay git** | `git init` + tag; documentar semillas y versiones |
| Método | Tres leyes de MAREA-2 | Logística y Student-t sobreviven a A2/A3; `sin_quorum` sólo en Fase A (5k/50k, vs FULL): 1/48/22, 1/60/10, 0/28/43, 1/27/43 | Evidencia preliminar a favor. Confirmar a 500k sobre `sin_colonias` antes de escribir el método |
| Presupuesto efectivo | MaxFES = 5k / 50k / 500k | **No es igual para todos.** Clásicos: (T+1)·N = **5040 / 50 080 / 500 080** (`competidores.py:48`, población inicial + T iteraciones). MAREA-2: N·T = **4960 / 49 920 / 499 920** (`experimento_config.py:157`). MAREA-2 corre con 160 evaluaciones menos: −3,2 % a 5k, despreciable a 500k | Sesgo conservador (contra la propuesta): se declara, no se corrige. L-SHADE/jSO reciben (T+1)·N como los clásicos. Chequeo barato opcional a 5k (Fase 2) |

---

## Fase 0 · Cerrar decisiones y congelar el algoritmo — HOY

- [ ] **A.** MAREA-2 sin_colonias como propuesta principal; FULL sólo en la tabla de ablación.
- [ ] **B.** Media como métrica principal (lo que pide el setup); mediana reportada como robustez. Declararlo en el texto.
- [ ] **C.** Reportar siempre rango de Friedman + W/T/L juntos; nombrar las familias de funciones donde gana DE en D10/500k.
- [ ] **D.** MAREA, MAREA_REL, MAREA-bio como línea base y ablación histórica ("lo que no funcionaba y por qué"), no como contendientes.
- [ ] **Tests pareados (cerrado):** rank-sum por función (51 vs 51 corridas, muestras independientes) para el W/T/L del §7; signed-rank pareado por problema sobre las 71 medias (Derrac et al.) como p-valor global por par de algoritmos y celda. Justificar en el texto que `seed = run` no empareja corridas entre algoritmos.
- [ ] **E.** Competidores nuevos: **L-SHADE y jSO** (ambos con resultados publicados en CEC2017 a 10 000·D, validables). EA4eig/AGSK sólo si sobra tiempo: portar desde MATLAB es más caro y más riesgoso.
- [ ] **Frase para el paper** justificando 11 competidores en vez de 20–25: cobertura de paradigmas (evolutivo DE/L-SHADE/jSO, swarm PSO/GWO/WOA/FA/ABC, inmune CLONALG), ausencia de implementaciones Python validadas de los ganadores CEC, prioridad a la fidelidad y reproducibilidad frente a la cantidad.
- [ ] **Presupuesto efectivo (cerrado):** los 7 clásicos consumen (T+1)·N = 5040 / 50 080 / 500 080; MAREA-2 consume N·T = 4960 / 49 920 / 499 920. La diferencia (160 evaluaciones) corre **en contra** de la propuesta, así que se declara como sesgo conservador y no se re-corre nada. L-SHADE y jSO reciben (T+1)·N, igual que los otros competidores. Los dos presupuestos van al setup y a la tabla de parámetros, y la nota se coloca junto al resultado de D50/5k (la única celda donde MAREA-2 no es 1º y donde un 3,2 % podría pesar).
- [ ] **`git init`** en el proyecto; commit del estado actual.
- [ ] **Añadir a `marea2.config()` el brazo `sin_colonias + sin_quorum`** (hoy no existe). Va **antes** del tag porque toca `marea2.py`.
- [ ] **Comprobar que `sin_colonias` no cambió — en el Mac mini:** una celda corta (p. ej. D10/5k, 5 corridas, mismas semillas) contra el JSON existente del protocolo; los valores deben ser idénticos bit a bit. Tiene que correr en el mini: en el MacBook sólo 2 de 20 corridas salen idénticas por diferencias de BLAS, y la prueba fallaría con el código correcto. Sin esto el tag no vale.
- [ ] Tag `marea2-sc-frozen` sobre `marea2.py` + `experimento_config.py`; guardar el hash SHA-256 de ambos en el archivo de decisiones.
- [ ] Congelar MAREA-2 SC: tabla de parámetros con la config única y lista de semillas. Desde aquí, cero cambios de diseño.
- [ ] **Lanzar en el mini la Fase 1a** (`sin_colonias + sin_quorum`): D10/D50 × 50k/500k, 15 corridas, mismas semillas que la A3.

**Listo cuando:** existe un archivo de decisiones (una página), el repo tiene tag con `sin_colonias` verificado idéntico en el mini, y la Fase 1a está corriendo.

---

## Fase 1a · Confirmar la tercera ley (quorum) — Mac mini, ~3–4 h de máquina

**Lo que ya hay (Fase A, `resultados_experimento/ablacion/`, 15 corridas):**

| `sin_quorum` vs FULL (W/T/L) | D10 5k | D10 50k | D50 5k | D50 50k |
|---|---|---|---|---|
| | 1/48/22 | 1/60/10 | 0/28/43 | 1/27/43 |

Quitar el quorum empeora en las 4 celdas, con mucha más fuerza en D50. Es evidencia preliminar a favor, pero con dos limitaciones que el paper no puede saltarse: (i) está medida a 5k/50k, y la memoria ya registró que el screening a 5k/50k no predijo el 500k (las colonias perjudicaban ahí y FULL fue 1º a 500k igualmente); (ii) está medida contra FULL, así que el efecto del quorum puede estar confundido con el de las colonias.

**Lo que falta:**

- [ ] `sin_colonias + sin_quorum` vs `sin_colonias`, D10/D50 × 50k/500k, 15 corridas, semillas de la A3 (pareable con la A3 y comparable con la fila de Fase A). Coste estimado: ~1 h D10/500k, ~2 h D50/500k, 50k despreciable.
- [ ] Rango de Friedman entre los brazos + W/T/L por rank-sum, como en la A3.
- [ ] Resultado esperado: se confirma. Entonces entra como tercera ley con su fila en la tabla de ablación (Fase A + esta confirmación).
- [ ] Si no se confirma a 500k sobre `sin_colonias`: el método se presenta con dos leyes y el trade-off μmax–biofilm como componente que ayuda a presupuesto corto y se agota a 500k, con ambos datos en la tabla. Es un resultado tan publicable como el positivo y refuerza la narrativa anti-metáfora.
- [ ] Actualizar `MAREA2_especificacion.md` con el resultado.

**Listo cuando:** la tabla de ablación tiene una fila por cada ley que el método afirma.

---

## Fase 1b · Implementar y validar L-SHADE y jSO — MacBook, 2–4 días, sin máquina de producción

- [ ] Revisar si `pyade` (tiene L-SHADE, jSO y variantes) sirve como base. Aunque se use, la validación contra CEC2017 sigue siendo obligatoria.
- [ ] **Interfaz:** los competidores nuevos reciben `maxfes` y derivan su propia población, no `(N, T)`. L-SHADE: N_init = 18·D con reducción lineal hasta 4. jSO: N_init = 25·log(D)·√D. Nada de `POP_SIZE = 80` aquí.
- [ ] Parámetros de los autores originales (memoria H = 5 en jSO, H = 6 en L-SHADE, archivo externo, p-best, etc.) sin re-tuning. Documentarlos en la tabla de parámetros.
- [ ] **Validación contra CEC2017 D50 a 500 000 exactos** (= 10 000·D, el presupuesto oficial; aquí se usa el nominal, no el efectivo del protocolo): ~25 corridas de cada uno sobre las 29 funciones de CEC2017 D50 y comparar media y mediana por función (los JSON ya guardan f − f*, así que la comparación es directa) con las tablas publicadas. Criterio: mismo orden de magnitud por función y mismo patrón de funciones resueltas (error < 1e-8). No se exige igualdad: las semillas son distintas. 51 corridas serían varias horas por algoritmo en el MacBook y no añaden nada a una validación.
- [ ] En el protocolo, pasarles **(T+1)·N = 5040 / 50 080 / 500 080**, el mismo presupuesto que los 7 clásicos, no el MaxFES nominal ni el N·T de MAREA-2.
- [ ] Chequeos de presupuesto en una corrida corta por celda: FES consumidos entre 0,98× y 1,02× de (T+1)·N.
- [ ] Verificar que la reducción de población recibe el MaxFES de **cada celda** (5k, 50k, 500k), no un valor fijo.

**Listo cuando:** ambos algoritmos reproducen el patrón de CEC2017 D50 y pasan el chequeo de presupuesto en las 6 celdas.

---

## Fase 1c · Protocolo de los competidores nuevos — Mac mini, ~2–4 días de máquina

Sólo el mini genera resultados (el MacBook da 2/20 corridas idénticas).

- [ ] Preflight: `verificar_entorno.py`, `correr_marea.py --check`, confirmar opfunu **1.0.4** (1.0.1 da no-finito en CEC2017 F17/F29 a D10).
- [ ] Mover cualquier JSON viejo fuera de `resultados_experimento/` (`run_experiment` reanuda y no recalcula).
- [ ] Orden de celdas: **D10/500k y D50/500k primero** (donde un DE moderno más puede cambiar el titular). Con dos algoritmos, las celdas de 500k de ambas dimensiones tardan más de un día; contar con primeros resultados en 1–2 días, no en horas.
- [ ] Luego 50k, luego 5k.
- [ ] **Validación por celda (algoritmo nuevo, no monotonicidad):**
  - presupuesto medido ≈ 1,0× de (T+1)·N (5040 / 50 080 / 500 080);
  - 3621 entradas (71 × 51);
  - ningún valor no finito.
- [ ] Guardar con el mismo esquema de nombres: `resultados_experimento/{ALG}_D{10,50}_FES{5000,50000,500000}.json`.

**Listo cuando:** 6 celdas × 51 corridas por competidor nuevo, las tres validaciones pasadas.

---

## Fase 2 · Análisis: pipeline conforme al setup — MacBook, en paralelo con 1a/1b, 3–5 días

Objetivo: un solo comando que lea todos los JSON y produzca todas las tablas y figuras, para que la Fase 3 sea sólo volver a correrlo.

- [ ] Lista de algoritmos leída de configuración, no fija en el código (hoy `generar_figuras.py`, `generar_pdfs_resultados.py` y `generar_pdfs_corridas.py` tienen los 11 fijos).
- [ ] Media como métrica principal; mediana como robustez.
- [ ] Regla de empate <1e-8 (CEC 2017).
- [ ] Rankings separados: por benchmark × D × MaxFES (CEC2014 D10/D50, CEC2017 D10/D50, **CEC2022 D10/D20** → 18 rankings) y agregados CEC_71_10 / CEC_71_50 × MaxFES (6). El conjunto "D50" se etiqueta en tablas y texto como "D50 (CEC2022 en D20)".
- [ ] Friedman + post-hoc de Shaffer, α = 0,05; % de competidores significativamente distintos por configuración.
- [ ] W/T/L por rank-sum contra **todos** los rivales (faltan PSO, CLONALG, MAREA, MAREABIO, MAREABIO_REL).
- [ ] Signed-rank pareado por problema (71 medias) por par de algoritmos y celda: p-valor global.
- [ ] SD de los rangos promedio por configuración.
- [ ] Varianza de la posición de cada algoritmo en los 6 rankings (§9, robustez).
- [ ] Tablas §8 en CSV y LaTeX: rangos por algoritmo, resumen de posiciones, Wilcoxon W/T/L, signed-rank global, % Friedman.
- [ ] Suplementario: media, mediana y SD de las 51 corridas por función × algoritmo × celda.
- [ ] Figuras regeneradas con MAREA-2 SC; FULL sólo en la figura/tabla de ablación. Reescribir título de fig2.
- [ ] Tabla de ablación de MAREA-2: colonias, Student-t vs gaussiana, logística vs Gompertz, **quorum (Fase A a 5k/50k + Fase 1a a 500k)**, memoria/lag/limit inertes, perturbación por eje y dimensión. Aclarar que `ruido_normal` no es gaussiano (el gaussiano real es `eje_gauss`).
- [ ] **Chequeo opcional del sesgo de presupuesto a 5k** (minutos de mini): re-correr MAREA-2 SC en D10/5k y D50/5k con T = 63 (5040 evaluaciones, mismas semillas, misma configuración) y comprobar que los puestos no cambian. Va al suplementario; no sustituye al protocolo. Cierra la objeción "su derrota en D50/5k es presupuesto" antes de que la haga un revisor.
- [ ] **Tabla de parámetros** de todos los algoritmos con fuente: declarar `POP_SIZE = 80` uniforme en los 7 clásicos como decisión de diseño (desviación del §4), gamma reescalado en FA (corrección de underflow de exp(−r²), desviación del §4), y población propia en L-SHADE/jSO.

**Control de regresión:** el pipeline debe reproducir la tabla 2.1 de la memoria (12 algoritmos, mediana) antes de cambiar nada.
**Listo cuando:** corre de una vez sobre los 13 algoritmos actuales y pasa el control de regresión.

---

## Fase 3 · Integrar y reevaluar — 1 día, al cerrar la Fase 1c

- [ ] Correr el pipeline con L-SHADE y jSO.
- [ ] Reescribir el titular con lo que salga. Si MAREA-2 deja de ser 1º, el paper pasa a "dónde compite y dónde no": es lo que el §9 pide de todos modos; la ablación y la narrativa anti-metáfora siguen intactas.
- [ ] Decisión C con datos: familias de funciones (unimodal / multimodal / híbrida / composición) donde gana DE, L-SHADE, jSO.
- [ ] Actualizar la memoria del proyecto.

---

## Fase 4 · Conformidad y reproducibilidad — 2 días

- [ ] Actualizar el doc de setup: 71 funciones y por qué; CEC2022 en D20 dentro del conjunto "D50"; los dos presupuestos efectivos (clásicos y L-SHADE/jSO: (T+1)·N = 5040 / 50 080 / 500 080; MAREA-2: N·T = 4960 / 49 920 / 499 920) con la nota de que el sesgo es contra la propuesta; 3 presupuestos con justificación (500k = convención CEC para 50D; 5M descartado por coste, 7 de 22 celdas); tests pareados y por qué; lista final de competidores y la frase que justifica 11.
- [ ] Tabla de parámetros y tamaños de población de todos los algoritmos (de la Fase 2), con las tres desviaciones del §4 declaradas.
- [ ] Paquete reproducible: repo con tag; semillas de las 51 corridas (documentar que protocolo y ablación usan semillas distintas: 0/1065 coincidencias); versiones exactas de implementaciones (opfunu 1.0.4, pyade u origen del código de L-SHADE/jSO, con su hash); JSON crudos.

---

## Fase 5 · Escritura — desde la semana 4

Orden sugerido (los datos mandan):

1. Resultados (§8) y ablación de MAREA-2.
2. Método: las tres leyes de Stehlik et al. 2025 (logística xmid 3,429 / s 0,784; Student-t ν = 3,05; trade-off μmax–biofilm, ec. 7 → regla 1/5 del paso, **con la evidencia de Fase A + Fase 1a**) y la regla de MAREA-2.
3. Sección "lo que no funcionaba y por qué": ec. (1) de AGEO infundada; σ sin escalar; colapso de fases de MAREA-bio; inversión bio→base.
4. Discusión (§9): dónde destaca y dónde no; componentes útiles por MaxFES; robustez de posición.
5. Introducción al final.
6. **Teorema de convergencia (decisión I):** reescribir `convergencia_proof_draft.md` para la regla de MAREA-2 SC, no la de AGEO. Va después de cerrar el texto del método; no bloquea nada anterior.

---

## Calendario

| Semana | Mini | MacBook |
|---|---|---|
| 1 | Fase 1a (`sin_colonias + sin_quorum`, ~3–4 h); luego libre | Fase 0 hoy (`git init`, brazo nuevo, verificación de `sin_colonias`, tag); Fase 1b (implementar/validar L-SHADE, jSO); Fase 2 |
| 2 | Fase 1c: protocolo L-SHADE y jSO, 500k primero | Fase 2 termina; análisis de 1a |
| 3 | Fase 1c termina | Fase 3; Fase 4 |
| 4+ | libre | Fase 5 |

El riesgo de calendario está en la Fase 1b (implementación y validación) y en la 1c (máquina). Todo lo demás se hace con los datos que ya existen.

---

## Reglas operativas (para no romper resultados)

- Sólo el Mac mini genera resultados; el MacBook implementa y analiza.
- opfunu 1.0.4, nunca 1.0.1. Correr `--check` y `verificar_entorno.py` antes de calcular.
- `early_stop` de MAREA desactivado por defecto; la validación por monotonicidad aplica **sólo** a re-corridas de MAREA con las mismas semillas, no a competidores nuevos.
- Competidores nuevos: validar presupuesto (~1,0× de (T+1)·N), 3621 entradas por celda y ausencia de no-finitos.
- Presupuestos efectivos: clásicos y competidores nuevos (T+1)·N = 5040 / 50 080 / 500 080; MAREA-2 N·T = 4960 / 49 920 / 499 920. No igualarlos re-corriendo: declararlos.
- Cualquier verificación de identidad bit a bit corre en el mini, nunca en el MacBook.
- Mover JSON viejos antes de relanzar (`run_experiment` reanuda).
- Resultados en `resultados_experimento/`; ablación en `resultados_experimento/ablacion/`; las descargas caen en la raíz: moverlas.
- Todo cambio de código pasa por commit; el tag `marea2-sc-frozen` no se mueve.
