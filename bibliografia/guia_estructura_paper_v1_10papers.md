# Cómo se escribe un paper que propone una metaheurística — síntesis de 10 papers

Fecha: 2026-10-04. PDFs en `bibliografia/metaheuristicas/` (fuera de git).

## Corpus

| # | Algoritmo | Revista, año | Inspiración |
|---|---|---|---|
| 01 | GWO | Advances in Engineering Software, 2014 | lobos |
| 02 | SMA | Future Generation Computer Systems, 2020 | moho mucilaginoso |
| 03 | HGS | Expert Systems with Applications, 2021 | hambre (etología) |
| 04 | RUN | Expert Systems with Applications, 2021 | Runge-Kutta ("beyond the metaphor") |
| 05 | INFO | Expert Systems with Applications, 2022 | media ponderada (sin metáfora) |
| 06 | HHO | Future Generation Computer Systems, 2019 | halcones |
| 07 | RIME | Neurocomputing, 2023 | escarcha (física) |
| 08 | HO | Scientific Reports, 2024 | hipopótamos |
| 09 | APO | Knowledge-Based Systems, 2024 | protozoos (Euglena) |
| 10 | SBOA | Artificial Intelligence Review, 2024 | ave secretaria |

Complementarios (los que leen los revisores críticos):
- Aranha et al. (2022), *Swarm Intelligence* 16:1–6, "the elephant in the room".
- Vermetten et al. (2024), arXiv:2402.09800, benchmarking a gran escala de heurísticas con metáfora.

Sesgo del corpus: 9 de 10 tienen como coautor a Mirjalili o a Heidari/Chen (GWO, HO y APO: Mirjalili; SMA, HGS, RUN, INFO, HHO y RIME: Heidari/Chen). Solo SBOA es independiente. Es representativo del género, pero ese grupo tiene un estilo propio.

---

## 1. El esqueleto común

Orden de secciones que se repite casi sin variación:

1. **Abstract.** Propuesta → mecanismos o fases → validación (nº de funciones y competidores) → resultado → aplicación → enlace al código.
2. **Introducción.** Límites de los métodos exactos o de gradiente → metaheurísticas → exploración/explotación → No Free Lunch → propuesta → contribuciones → párrafo de organización.
3. **Revisión de literatura.** Taxonomía por fuente de inspiración. En 5/10 es una sección propia; en los otros 5 va dentro de la introducción.
4. **Inspiración.** 1–2 páginas con figuras, fotos o esquemas.
5. **Modelo matemático.** Inicialización, fases de exploración y explotación, pseudocódigo, diagrama de flujo y complejidad O(·).
6. **Experimentos.**
   - Setup.
   - Análisis cualitativo.
   - Benchmarks, uno por suite.
   - Estadística.
   - Convergencia.
   - Sensibilidad de parámetros.
   - (Escalabilidad y tiempo, a veces.)
7. **Problemas de ingeniería con restricciones**, a veces con una aplicación real.
8. **Conclusiones y trabajo futuro.**

## 2. Puntos comunes (frecuencia en los 10)

| Elemento | Frecuencia | Observaciones |
|---|---|---|
| Párrafo "The rest of the paper is organized as follows" | 10/10 | En 3 no coincide con las secciones reales |
| Pseudocódigo | 10/10 | |
| Fases explícitas de exploración y explotación | 10/10 | Casi siempre asignadas por narrativa, no medidas |
| Problemas de ingeniería | 10/10 | Recipiente a presión, viga soldada y resorte en casi todos |
| 30 corridas por función | 10/10 | Solo 30, ninguno usa 51 (CEC) |
| Teorema No Free Lunch como motivación | 9/10 | Solo RUN no lo usa |
| Complejidad computacional O(·) | 9/10 | En 4 la fórmula está mal |
| Código publicado o prometido | 9/10 | Datos crudos: 0/10 |
| Taxonomía de familias (evolutivo, enjambre, física, humano) | 10/10 | Tabla de ventajas y desventajas solo en SBOA |
| Curvas de convergencia | 9/10 | |
| Análisis cualitativo de 4 paneles: historial, trayectoria, fitness medio, convergencia | 8/10 | Plantilla de Mirjalili |
| Diagrama de flujo | 8/10 | |
| Wilcoxon | 8/10 | Rank-sum o signed-rank, casi nunca corregido |
| Friedman | 8/10 | |
| Sección de inspiración con figuras | 8/10 | INFO y RUN la sustituyen por una base matemática |
| Sensibilidad de parámetros | 7/10 | Siempre ajustada sobre las mismas funciones de evaluación |
| Lista explícita de contribuciones en la introducción | 4/10 | RIME, HO, APO y SBOA |
| Tiempo de cómputo | 5/10 | |
| Escalabilidad en dimensión | 4/10 | |
| Post-hoc de Friedman | 3/10 | Holm (SMA), Bonferroni-Dunn (INFO), Nemenyi (HO) |
| Boxplots | 2/10 | |
| Porcentajes de exploración y explotación (diversidad) | 2/10 | APO y SBOA |
| Comparación con L-SHADE o variantes | 2/10 | HGS y SBOA; jSO en 0/10 |
| Presupuesto igualado en FEs | 3/10 | Los otros 7 comparan por iteraciones |
| **Ablación por componente** | **0/10** | |
| **Datos biológicos medidos dentro de las ecuaciones** | **0/10** | APO cita datos de Euglena que no usa |
| Sección honesta de limitaciones | ~0/10 | HO tiene una genérica |

## 3. Lo que hay que copiar

- **Abstract en ese orden**, cerrado con el enlace al código o los datos.
- **Contribuciones en viñetas** al final de la introducción (como RIME, HO y APO), más el párrafo de organización.
- **Taxonomía breve** de metaheurísticas para situar la propuesta. Los críticos piden situarla también por mecanismo, no solo por inspiración (ver §5).
- **Pseudocódigo, diagrama de flujo y complejidad O(·) correcta.** Que esté bien es un diferenciador barato: 4 de 9 la tienen mal.
- **Análisis cualitativo** (historial de búsqueda, trayectoria, convergencia). Los revisores ya lo esperan.
- **Estadística por capas:** Friedman → post-hoc → Wilcoxon por pares, con conteo +/=/−. El diagrama CD (HO) es muy legible.
- **Sensibilidad de parámetros**, pero sin repetir el defecto general de ajustarla sobre las mismas funciones.
- **Escalabilidad** (D10/D50) y **tiempo de cómputo**.
- **Código público** (9/10 lo tienen; no tenerlo sería raro).

## 4. Lo que NO hay que copiar (errores recurrentes)

1. **Comparar por iteraciones** con algoritmos que gastan distinto nº de FEs por iteración (7/10). HHO, RUN y SBOA gastan 2–3 FEs por agente e iteración.
2. **Presupuestos muy por debajo del estándar CEC (10000·D):** 30k FEs en CEC2014 D30, 60k en D100.
3. **Competidores débiles o variantes de la misma familia**, sin L-SHADE, jSO ni CMA-ES bien configurado. Cuando HGS incluye L-SHADE le fija N=30, lo que anula su reducción de población.
4. **Friedman sin post-hoc** y Wilcoxon sin corrección por comparaciones múltiples.
5. **Funciones con el óptimo en el origen** (F1–F13 clásicas) sin controlar el sesgo al centro. Varios algoritmos (GWO, SMA, HGS, HO) contraen hacia 0.
6. **Ingeniería con soluciones infactibles.** Los agentes lo verificaron en 6/10: variables discretas reportadas como continuas, restricciones violadas o costos por debajo del óptimo conocido. Rivales copiados de la literatura con presupuestos distintos y solo el mejor valor.
7. **Metáfora que no se modela o se abandona.** HO lo admite ("would reduce the performance"); en GWO los roles ω y δ no aparecen en las ecuaciones.
8. **Operadores prestados sin decirlo.** Mutación diferencial DE (SMA, RUN, INFO, APO), Lévy y browniano (SBOA), selección greedy de DE (RIME).
9. **Inconsistencias internas:** nº de competidores y de problemas distinto entre abstract y texto, columnas imposibles (SD negativa en GWO, MED > media en SMA), erratas.
10. **Conclusiones triunfales sin limitaciones.**

## 5. Lo que exigen los revisores críticos

### Aranha et al. (2022)

Una revista no debería publicar el algoritmo salvo que:
1. esté descrito en **terminología estándar de optimización**;
2. aporte **conceptos útiles y nuevos** al campo;
3. la metáfora esté motivada **sobre una base científica**;
4. la comparación sea **justa, contra el estado del arte** y con buenas prácticas de benchmarking.

Rechazan que metáfora, modelo e implementación sean "three (almost completely) different things". También rechazan enviar el trabajo a revistas del área de la metáfora (p. ej. microbiología) para esquivar a los revisores de optimización.

### Vermetten et al. (2024)

- 294 implementaciones con metáfora en BBOB. Muchas pierden contra RandomSearch.
- Arriba quedan las referencias (BIPOP-CMA-ES, DE, L-SHADE).
- Los metafóricos destacan sobre todo **con presupuesto bajo**.

Recomiendan:
- incluir **RandomSearch** como línea base;
- publicar el código;
- presentar el algoritmo como modificación de uno existente y **medir cada cambio** (ablación);
- declarar la paridad de ajuste entre algoritmos;
- usar más de una suite;
- reportar análisis anytime además del de presupuesto fijo;
- valorar el perfil de especialista.

---

## 6. Dónde queda MAREA-2 frente a esto

### Ya lo tenemos y casi nadie en el corpus lo tiene (son nuestros argumentos)

| Diferenciador | Corpus | MAREA-2 |
|---|---|---|
| Leyes biológicas **medidas** dentro de las ecuaciones | 0/10 | Logística (xmid 3,429; s 0,784), Student-t ν=3,05, trade-off µmax–biofilm convertido en regla 1/5 |
| Ablación por componente | 0/10 | `t11_ablacion`: 22 comparaciones |
| L-SHADE **y** jSO como rivales | 2/10 y 0/10 | Sí, con sus parámetros publicados y presupuesto igualado |
| Presupuesto en FEs al nivel CEC (10000·D en D50/500k) | 3/10 por FEs | 5k / 50k / 500k, déficit de 80 FEs declarado |
| 51 corridas (protocolo CEC) | 0/10 (todos 30) | 51 |
| Tres suites: CEC2014, CEC2017, CEC2022 | Máximo 2 | 3 suites, 71 funciones |
| Friedman + Shaffer + Nemenyi + Wilcoxon con Holm | 3/10 tienen post-hoc | `t4`–`t9`, `fig_cd_nemenyi` |
| Reportar dónde pierde | ~0/10 | Pierde con L-SHADE/jSO a 50k/500k: se dice |
| Datos crudos públicos | 0/10 | JSON de 3621 corridas por celda |

### Lo que el género espera y aún no tenemos

| Elemento | Estado | Costo |
|---|---|---|
| Curvas de convergencia | No guardamos historia (`record_history=False`) | Corrida nueva en el mini, con historia y pocas funciones |
| Análisis cualitativo (historial y trayectoria en 2D) | No hecho | Bajo: se puede correr en el MacBook (es solo ilustración, no resultado) |
| Boxplots | No hecho | Bajo: sale de los JSON |
| Complejidad O(·) de MAREA-2 | No escrita | Bajo: análisis del código |
| Tiempo de cómputo | No medido | Bajo–medio: en el mini |
| Sensibilidad de parámetros | Solo la ablación | Medio: barrido en el mini, idealmente en otra suite que la de evaluación |
| RandomSearch como línea base | No | Bajo en el mini (lo pide Vermetten) |
| Problemas de ingeniería con restricciones | No | Medio: 10/10 los traen; depende de la revista (ver §7) |
| Diversidad / % exploración-explotación | No | Medio: necesita historia |
| Código en Zenodo con DOI | Repo privado | Bajo, al enviar |

## 7. Esqueleto propuesto para el paper de MAREA-2

1. **Abstract**, en el orden del §1. Cierra con el código y los datos.
2. **Introduction**
   - Problema: optimización de caja negra; proliferación de metáforas sin aporte (Sörensen 2015; Aranha 2022). Usar el NFL como contexto, no como licencia.
   - Gap: los algoritmos bioinspirados describen comportamiento, pero ninguno traduce **leyes cuantitativas medidas** del organismo a reglas de búsqueda ni las somete a ablación.
   - Contribuciones en viñetas:
     1. tres leyes medidas de Cobetia → tres reglas estándar;
     2. ablación;
     3. evaluación con 71 funciones, 51 corridas, 3 presupuestos y estado del arte;
     4. el perfil "dónde compite y dónde no".
   - Párrafo de organización.
3. **Related work**, con dos ejes:
   - la taxonomía por inspiración, breve;
   - **los mecanismos equivalentes**: reducción de población (LPSR de L-SHADE), mutación de cola pesada (Cauchy y Lévy en la programación evolutiva rápida), regla 1/5 (Rechenberg) y DE de estado del arte.

   Situar cada pieza de MAREA-2 frente a su equivalente estándar (exigencia 1 de Aranha).
4. **De la biología a las reglas**
   - Cobetia y el paper de Frontiers.
   - Tabla de trazabilidad en tres columnas: **ley medida (fuente) → ecuación → componente del código**. Responde directamente a "three different things".
   - Sin vocabulario metafórico operativo.
5. **MAREA-2**
   - Algoritmo en terminología estándar, pseudocódigo, diagrama de flujo, complejidad O(·) y parámetros (`t12`).
   - Transferencia horizontal parcial: qué es en términos estándar (cruce por coordenadas con fracción adaptativa) y qué la distingue del CR adaptativo de SHADE.
6. **Diseño experimental**
   - Suites, D, presupuestos (5k/50k/500k en FEs), 51 corridas, rivales con sus parámetros publicados, déficit de 80 FEs.
   - Métricas: media como principal, mediana como robustez, regla de 1e-8.
   - Tests: Friedman, Shaffer, Nemenyi y Wilcoxon con Holm.
   - Separar el conjunto de diseño (ablación con 15 semillas en CEC2014) del de evaluación.
7. **Resultados**
   - Rangos (`t1`, `fig_rangos`), por benchmark (`t2`), CD (`fig_cd_nemenyi`), W/T/L (`t7`, `fig_wtl`), por clase de función (`t8`) y robustez (`t3`, `t10`).
   - Convergencia y boxplots, cuando los tengamos.
8. **Ablación** (`t11`): qué aporta cada ley y qué es inerte.
9. **Discusión**
   - Perfil de especialista: gana a 5k y en D50 5k a todos; detrás de L-SHADE/jSO a 50k/500k.
   - Familias de funciones donde pierde: híbridas y unimodales.
   - Limitaciones explícitas.
10. **Conclusiones y trabajo futuro.**
11. **Disponibilidad de código y datos**: Zenodo.

## 8. Qué no afirmar

- "Supera al estado del arte" o "nuevo paradigma".
- "La bacteria optimiza" o "resuelve problemas".
- Términos metafóricos como nombres de operadores sin su equivalente estándar al lado.
- El resultado a 5k presentado como superioridad general.

## 9. Decisiones pendientes (del usuario y su profesor)

1. **Revista objetivo.** Swarm and Evolutionary Computation, IEEE TEVC, Evolutionary Computation, Applied Soft Computing… Decide si hacen falta problemas de ingeniería: todas las del corpus los traen, pero las revistas de computación evolutiva pura los exigen menos.
2. **Corridas nuevas en el mini:** curvas de convergencia, tiempo, RandomSearch y sensibilidad. ¿Cuáles valen el costo?
3. **Ingeniería sí o no.** Si sí: 3–5 problemas clásicos con manejo de restricciones declarado, estadística completa y los rivales corridos por nosotros, no copiados de la literatura.
