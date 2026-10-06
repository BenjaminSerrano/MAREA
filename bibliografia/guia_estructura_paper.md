# Cómo se escribe un paper que propone una metaheurística — síntesis de 20 papers

Versión 2, 2026-10-06. Reemplaza la versión de 10 papers del 2026-10-04 (respaldo en `guia_estructura_paper_v1_10papers.md`).
PDFs en `bibliografia/metaheuristicas/` (fuera de git). Una ficha JSON por paper en `bibliografia/fichas/` (esquema en `fichas/FICHA.md`); todas las cifras de esta guía salen de esas fichas.

## Corpus

La versión 1 tenía 9 de 10 papers con Mirjalili o Heidari/Chen como coautores. Se dejaron 5 de ellos, se sacaron 5 (SMA, HGS, INFO, RIME, HO; en `metaheuristicas/descartados/`) y se agregaron 15 de grupos distintos, con un máximo de un paper por grupo. Ahora Mirjalili aparece en 3 de 20 y Heidari/Chen en 2 de 20.

| # | Algoritmo | Autores (primero) | Revista, año | Tipo | Inspiración |
|---|---|---|---|---|---|
| 01 | GWO | Mirjalili | Adv. Eng. Software, 2014 | metáfora | lobos |
| 04 | RUN | Ahmadianfar | ESWA, 2021 | sin metáfora* | Runge-Kutta |
| 06 | HHO | Heidari | FGCS, 2019 | metáfora | halcones |
| 09 | APO | Wang | KBS, 2024 | microbiano | protozoo *Euglena* |
| 10 | SBOA | Fu | AIR, 2024 | metáfora | ave secretaria |
| 11 | JADE | Zhang | IEEE TEVC, 2009 | sin metáfora | — |
| 12 | L-SHADE | Tanabe | IEEE CEC, 2014 | sin metáfora | — |
| 13 | jSO | Brest | IEEE CEC, 2017 | sin metáfora | — |
| 14 | CMA-ES | Hansen | Evol. Comput., 2001 | sin metáfora | — |
| 15 | BFO | Passino | IEEE Control Syst. Mag., 2002 | microbiano | quimiotaxis de *E. coli* |
| 16 | BCO | Niu | Discrete Dyn. Nat. Soc., 2012 | microbiano | ciclo de vida de *E. coli* |
| 17 | COA | Dehghani | KBS, 2023 | metáfora | coatíes |
| 18 | GJO | Chopra | ESWA, 2022 | metáfora | chacal dorado |
| 19 | KOA | Abdel-Basset | KBS, 2023 | metáfora | leyes de Kepler (física) |
| 20 | ChOA | Khishe | ESWA, 2020 | metáfora | chimpancés |
| 21 | RTH | Ferahtia | Sci. Rep., 2023 | metáfora | halcón de cola roja |
| 22 | GTO | Sadeeq | IEEE Access, 2022 | metáfora | pez jurel gigante |
| 23 | EHO | Al-Betar | AIR, 2024 | metáfora | manadas de wapitíes |
| 24 | BKA | Wang | AIR, 2024 | metáfora | milano de alas negras |
| 25 | EVO | Azizi | Sci. Rep., 2023 | metáfora | física nuclear |

\* RUN no usa metáfora, pero sigue la plantilla del grupo Heidari (NFL, taxonomía, fases, ingeniería). En las tablas cuenta como "sin metáfora"; cuando hablo del **género** me refiero a los 16 papers con metáfora, microbianos y RUN, y cuando hablo de **computación evolutiva (CE)** me refiero a JADE, L-SHADE, jSO y CMA-ES.

Grupos: metáfora 12, microbianos 3, sin metáfora 5.

Complementarios (los que leen los revisores críticos):
- Aranha et al. (2022), *Swarm Intelligence* 16:1–6, "the elephant in the room".
- Vermetten et al. (2024), arXiv:2402.09800, benchmarking a gran escala de heurísticas con metáfora.

Sesgos que quedan:
- **Acceso abierto:** de los 15 nuevos, 11 salieron de acceso abierto o de la página del autor y 4 (COA, GJO, KOA, ChOA) los bajó el usuario con acceso UV. Los papers cerrados de Elsevier quedan subrepresentados.
- **Versiones:** GWO, RUN y HHO son preprints; L-SHADE, jSO, BFO y CMA-ES son versiones de autor.
- **Muestra, no censo:** las frecuencias describen este corpus, no el género entero.

---

## 1. Dos maneras de escribir el paper

El hallazgo principal de la ampliación es que el corpus se parte en dos estilos que casi no se mezclan.

| | Género (16) | Computación evolutiva (4) |
|---|---|---|
| Motivación | NFL (12/16) y "las metaheurísticas son útiles" | Un defecto concreto del algoritmo padre |
| Justificación del diseño | Historia natural → ecuación | Razonamiento por componente, requisitos, comparación con el padre |
| Mecanismos heredados | 13/16 dejan al menos uno sin declarar | 4/4 declaran y citan cada pieza |
| Ablación | 0/16 | 3/4 (JADE, L-SHADE, CMA-ES) |
| Rivales | Muchos (mediana 10, hasta 32), sobre todo otros metafóricos | Pocos (2–8), el padre y el estado del arte |
| Presupuesto | Iteraciones (9/15 con dato) | FEs (3/4; JADE usa generaciones fijas) |
| Ingeniería | 14/16 | 0/4 |
| Extensión (págs. del PDF) | mediana ~30, hasta 102 | 8–37 |
| Tono | Promocional ("outstanding", "beyond a reasonable doubt") | Técnico; CMA-ES opina en primera persona |

MAREA-2 tiene que estar en medio: se inspira en Cobetia (género) pero quiere ser leída por revisores de CE. La propuesta del §7 toma la forma del género (secciones que el lector espera) y la argumentación de CE (requisitos, ablación, mecanismos declarados).

## 2. El esqueleto común del género

Orden de secciones que se repite casi sin variación en los 16:

1. **Abstract.** Propuesta → inspiración → fases → validación (nº de funciones y rivales) → resultado → ingeniería → (código).
2. **Introducción.** Límites de los métodos exactos → metaheurísticas → exploración/explotación → NFL → propuesta → contribuciones → organización.
3. **Revisión de literatura** con taxonomía por fuente de inspiración (sección propia en 10/20).
4. **Inspiración**, 1–3 páginas con figuras (16/16 de los que tienen metáfora o microbio).
5. **Modelo matemático:** inicialización, fases de exploración y explotación, pseudocódigo, diagrama de flujo, O(·).
6. **Experimentos:** setup, análisis cualitativo, benchmarks, estadística, convergencia, (sensibilidad, escalabilidad, tiempo).
7. **Ingeniería con restricciones**, a veces con una aplicación real.
8. **Conclusiones y trabajo futuro.**

Los papers de CE son más cortos y van directo: contexto del padre → cambio propuesto → por qué → experimentos → (limitaciones).

## 3. Frecuencias

`n/N`: N excluye los papers donde el campo no aplica o no se pudo determinar. Columnas: M = metáfora (12), B = microbianos (3), S = sin metáfora (5).

### Estructura y presentación

| Elemento | M | B | S | Total | Observaciones |
|---|---|---|---|---|---|
| Párrafo de organización | 12/12 | 3/3 | 4/5 | 19/20 | No coincide con las secciones reales en 4 (RUN, HHO, COA, RTH) |
| Pseudocódigo | 11/12 | 3/3 | 4/5 | 18/20 | |
| Sección de inspiración con figuras | 12/12 | 3/3 | 1/1 | 16/16 | |
| Taxonomía de familias | 12/12 | 1/3 | 1/5 | 14/20 | Tabla de ventajas y desventajas solo en SBOA |
| NFL como motivación | 11/12 | 1/3 | 0/5 | 12/20 | Ningún paper de CE lo usa |
| Fases explícitas de exploración y explotación | 10/12 | 1/3 | 1/5 | 12/20 | Medidas solo en 2/20 (APO, SBOA) |
| Revisión de literatura como sección propia | 7/12 | 0/3 | 3/5 | 10/20 | |
| Complejidad O(·) | 7/12 | 1/3 | 2/5 | 10/20 | Correcta en 5; mal en 2 (RUN, EVO); dudosa en 3 (HHO, KOA, BKA) |
| Lista de contribuciones | 6/12 | 2/3 | 1/5 | 9/20 | |
| Diagrama de flujo | 6/12 | 2/3 | 1/5 | 9/20 | |
| Metáfora modelada por completo en las ecuaciones | 0/12 | 0/3 | — | 0/15 | Todas "parcial": partes de la historia sin ecuación |
| **Datos biológicos medidos dentro de las ecuaciones** | **0/12** | **0/3** | — | **0/15** | BFO y APO describen datos medidos, pero no los usan en las reglas |

### Diseño experimental

| Elemento | M | B | S | Total | Observaciones |
|---|---|---|---|---|---|
| Compara por FEs | 4/12 | 1/2 | 3/5 | 8/19 | APO, COA, KOA, EVO + L-SHADE, jSO, CMA-ES. 10 comparan por iteraciones, 2 ambiguos |
| Presupuesto al nivel CEC (10000·D) | 3/11 | 1/3 | 2/4 | 6/18 | APO, COA, KOA, EHO, L-SHADE, jSO |
| Corridas | | | | | 30 en 11/20; 51 en 3 (L-SHADE, jSO, COA); 50 en 2; 100, 20, 10 y 1 en uno cada uno |
| Parámetros de los rivales publicados | 9/12 | 1/2 | 3/5 | 13/19 | |
| Rivales: L-SHADE canónico | | | | 1/19 | Solo jSO; SBOA y KOA usan variantes (cnEpSin, SPACMA) |
| Rivales: jSO | | | | 0/19 | |
| Rivales: CMA-ES | | | | 4/19 | GWO, RUN, L-SHADE, EHO |
| Rivales: RandomSearch | | | | 0/20 | |
| Escalabilidad en dimensión | 7/12 | 0/3 | 4/5 | 11/20 | |
| Tiempo de cómputo | 4/12 | 1/3 | 3/5 | 8/20 | EVO usa el procedimiento T0/T1/T2 de CEC |
| Sensibilidad de parámetros | 3/12 | 1/3 | 2/5 | 6/20 | Las 6 sobre las mismas funciones de evaluación |
| Ajuste en un conjunto distinto del de evaluación | | | | 1/20 | Solo L-SHADE (ParamILS en CEC2013) |
| **Ablación por componente** | **0/12** | **0/3** | **3/5** | **3/20** | JADE, L-SHADE, CMA-ES |

### Estadística y figuras

| Elemento | M | B | S | Total | Observaciones |
|---|---|---|---|---|---|
| Wilcoxon (rank-sum o signed-rank) | 9/12 | 1/3 | 2/5 | 12/20 | En 5 el test está mal rotulado o mal leído (ver §5) |
| Friedman | 3/12 | 1/3 | 1/5 | 5/20 | |
| Sin ningún test estadístico | | | | 6/20 | GWO, JADE, CMA-ES, BFO, BCO; RTH declara ANOVA pero no lo muestra |
| **Post-hoc de Friedman** | | | | **1/20** | Solo EHO (Holm con algoritmo de control) |
| Corrección por comparaciones múltiples | | | | 1/20 | EHO |
| Conteo W/T/L | 2/12 | 1/3 | 2/5 | 5/20 | |
| Curvas de convergencia | 12/12 | 2/3 | 3/5 | 17/20 | **Eje en FEs solo en KOA (1/17)** |
| Análisis cualitativo (historial, trayectoria, fitness medio) | 6/12 | 2/3 | 1/5 | 9/20 | |
| Boxplots | 6/12 | 0/3 | 1/5 | 7/20 | |
| Diversidad o % de exploración | 1/12 | 1/3 | 0/5 | 2/20 | APO, SBOA |

### Ingeniería, reproducibilidad y honestidad

| Elemento | M | B | S | Total | Observaciones |
|---|---|---|---|---|---|
| Problemas de ingeniería con restricciones | 12/12 | 1/3 | 1/5 | 14/20 | Recipiente a presión, viga soldada y resorte en casi todos |
| Aplicación real además de los clásicos | 5/12 | 1/3 | 0/5 | 6/20 | GWO, APO, GJO, KOA, ChOA, RTH |
| Rivales de ingeniería copiados de la literatura | 5/11 | 0/1 | 1/1 | 6/13 | |
| **Algún problema verificable en la ingeniería** | | | | **13/14** | Solo APO sale limpio (ver §5) |
| Código público | 4/12 | 2/3 | 2/5 | 8/20 | |
| **Datos crudos públicos** | | | | **0/18** | |
| Dice en el texto dónde pierde | 8/12 | 2/2 | 5/5 | 15/19 | Casi siempre de pasada |
| **Sección de limitaciones honesta** | 0/12 | 1/3 | 2/5 | **3/20** | L-SHADE, CMA-ES, BFO. Genérica en 6; ausente en 11 |

## 4. Cómo justifican el diseño los papers sin metáfora

Es lo que más le sirve a MAREA-2, porque es lo que el género no hace y lo que piden Aranha y Vermetten.

- **JADE: razonamiento por componente y ablación.** Cada decisión trae su porqué: la media de Lehmer corrige el sesgo hacia F pequeños; Cauchy en F diversifica; la vida media 1/c fija la memoria. Después quita un componente a la vez (rand-JADE, nona-JADE, con y sin archivo) y explica cómo interactúan.
- **L-SHADE: comparación con el padre y ajuste separado.** Aísla lo nuevo (LPSR) comparando contra SHADE 1.1 sin él. Ajusta parámetros con ParamILS en CEC2013 F21–F25, que no son las funciones de evaluación (CEC2014), y publica rangos, valores por defecto y afinados. Tiene una sección de limitaciones.
- **jSO: diferencias marcadas.** Marca con ⇐ en el pseudocódigo exactamente qué líneas cambian respecto a iL-SHADE y lista los parámetros que difieren. Justifica los valores de forma empírica y lo dice ("without fine tuning").
- **CMA-ES: requisitos primero, verificación después.** Enuncia los requisitos del diseño antes del algoritmo (adaptación, rendimiento, invariancias, estacionariedad bajo selección aleatoria) y los comprueba uno por uno. PATH-ES (CMA-ES sin adaptar la covarianza) funciona como ablación. Mide FEs hasta un objetivo frente a n = 5…320. Tiene limitaciones y dice dónde no gana.

Dos papers del género aportan piezas útiles para lo mismo:
- **BFO:** lista explícitamente qué simplifica el modelo frente a la biología real (largos de carrera exponenciales, memoria de 4 s, asincronía, consumo de nutrientes) y da guías razonadas del efecto de cada parámetro.
- **RTH (Tabla 1) y EHO:** traducen la metáfora a conceptos genéricos (elitismo, selección, autoadaptación; líderes/seguidores/hijos). EHO además reconoce que su operador equivale a (μ+λ), ruleta y coeficientes tipo PSO.

## 5. Errores recurrentes (no copiar)

Cada ficha cita página. Lista no exhaustiva.

1. **Comparar por iteraciones con algoritmos que gastan distinto nº de FEs.** 10 de 19 comparan por iteraciones. Gastan más de 1 FE por agente e iteración sin declararlo: RUN (1–3), HHO (hasta 2), SBOA (2), GTO (3); en RTH y BKA no se puede determinar. APO hace lo correcto: ajusta la población de los rivales que gastan 2 o 4 FEs por individuo.
2. **Presupuestos bajo el estándar CEC:** 12 de 18 están por debajo de 10000·D (p. ej. GJO usa 6000 FEs).
3. **Rivales débiles o sin el estado del arte.** L-SHADE canónico en 1/19, jSO en 0/19, RandomSearch en 0/20.
4. **Estadística sin post-hoc ni corrección** (19/20 sin post-hoc), y **tests mal rotulados o mal leídos**: COA llama al rank-sum "test de medias"; en GTO los p-valores mínimos son los de un signed-rank pero lo llama rank-sum; en BKA es al revés; EVO invierte la lógica de normalidad de Kolmogorov-Smirnov; EHO lee p = 1.27e-10 como "mayor que 0,05".
5. **Estadísticas imposibles** en 6/20: desviaciones negativas (GWO), media en el óptimo con desviación no nula (HHO), media menor que el mínimo (COA, EHO), desviación muy superior al rango (ChOA, BCO).
6. **Conteos que no cuadran** entre abstract y texto en 5/20: GWO (rivales), SBOA (15/14/11), COA (11 vs 19 rivales), KOA (8 vs 9 problemas; 23 vs 12 funciones), EHO (9/10/12 rivales).
7. **Ingeniería que no resiste un recálculo** (13 de 14):
   - Variables enteras o discretas reportadas como continuas: GWO, HHO, SBOA, GJO, GTO, EVO, EHO.
   - Restricción violada: COA (recipiente, g3 = +427,8), HHO (embrague bajo t_min).
   - El costo reportado no sale de las variables publicadas: KOA (ninguna fila del tren de engranajes), GTO (voladizo), RUN (tres costos para un diseño), EHO (variables ×16).
   - Rivales en otra formulación (continua vs discreta del recipiente): GWO, HHO, GJO, EVO.
   - Sin vectores de diseño, inverificable: ChOA, RTH.
   - "Nuevo óptimo" falso: BKA (5887,37 > 5885,33), GTO.
8. **Mecanismos prestados sin decirlo** (13/16 del género): GWO (atracción tipo PSO), COA (factor de enseñanza de TLBO), GJO (energía de escape de HHO), KOA (masas de GSA), RTH (estructura de Bald Eagle Search, sin citarlo), BKA (presenta la mutación de Cauchy como propia), BCO (PSO), EVO (cruce con el mejor). SBOA y BFO sí los declaran.
9. **Metáfora que no se modela:** 0/15 la modelan por completo. ChOA trae ADN y filogenia que no entran en ninguna ecuación; en EHO la ruleta da más seguidores a los peores líderes, al revés de la historia.
10. **Ajustar o analizar la sensibilidad sobre las mismas funciones de evaluación** (6/6 de los que lo hacen).
11. **Convergencia por iteraciones** (16/17), que no es comparable entre algoritmos con distinto gasto por iteración.
12. **Conclusiones triunfales sin limitaciones** (ninguno de los 12 con metáfora tiene una sección honesta).

## 6. Lo que exigen los revisores críticos

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

## 7. Dónde queda MAREA-2 frente a esto

### Ya lo tenemos y casi nadie en el corpus lo tiene (son nuestros argumentos)

| Diferenciador | Corpus | MAREA-2 |
|---|---|---|
| Leyes biológicas **medidas** dentro de las ecuaciones | 0/15 | Logística (xmid 3,429; s 0,784), Student-t ν=3,05, trade-off µmax–biofilm convertido en regla 1/5 |
| Ablación por componente | 0/16 en el género; 3/4 en CE | `t11_ablacion`: 22 comparaciones |
| L-SHADE **y** jSO como rivales | 1/19 y 0/19 | Ambos, con sus parámetros publicados y presupuesto igualado |
| Post-hoc de Friedman con corrección | 1/20 | Shaffer + Nemenyi + signed-rank con Holm (`t4`–`t9`, `fig_cd_nemenyi`) |
| Presupuesto en FEs al nivel CEC | 6/18 | 5k / 50k / 500k; a 500k en D50 es 10000·D. Déficit de 80 FEs declarado |
| 51 corridas (protocolo CEC) | 3/20 | 51 |
| Tres presupuestos por dimensión | ~0/20 (KOA usa dos de forma incidental) | 5k / 50k / 500k en D10 y D50: el perfil "fuerte a presupuesto bajo" que describe Vermetten |
| Sección de limitaciones honesta | 3/20 (ninguno con metáfora) | Pierde con L-SHADE/jSO a 50k/500k: se dice con números |
| Datos crudos públicos | 0/18 | JSON de 3621 corridas por celda |

Ya no es diferenciador:
- **Número de suites.** KOA usa cuatro (CEC2014/17/20/22, todas en D10). Las nuestras son 3 suites × 2 dimensiones × 3 presupuestos; el argumento es la combinación, no el número.
- **"Decir dónde pierde."** 15/19 lo mencionan de pasada. Lo raro es una sección de limitaciones con números.

### Lo que el género espera y aún no tenemos

| Elemento | Corpus | Estado | Costo |
|---|---|---|---|
| Curvas de convergencia | 17/20 (en FEs: 1/17) | No guardamos historia (`record_history=False`) | Corrida nueva en el mini, en FEs: sería diferenciador |
| Análisis cualitativo (historial y trayectoria en 2D) | 9/20 | No hecho | Bajo: MacBook (es ilustración, no resultado) |
| Boxplots | 7/20 | No hecho | Bajo: sale de los JSON |
| Complejidad O(·) correcta | 5/20 correcta | No escrita | Bajo: análisis del código |
| Tiempo de cómputo | 8/20 | No medido | Bajo–medio: en el mini, idealmente con T0/T1/T2 de CEC (como EVO) |
| Sensibilidad de parámetros | 6/20, todas en las mismas funciones | Solo la ablación | Medio: barrido en el mini en una suite distinta (como L-SHADE) |
| RandomSearch como línea base | 0/20 | No | Bajo en el mini |
| Problemas de ingeniería con restricciones | 14/16 del género, 0/4 de CE | No | Medio: depende de la revista. Si se hacen, verificar factibilidad ya es diferenciador (13/14 fallan) |
| Diversidad / % exploración | 2/20 | No | Medio: necesita historia |
| Código en Zenodo con DOI | 8/20 código | Repo privado | Bajo, al enviar |

## 8. Esqueleto propuesto para el paper de MAREA-2

1. **Abstract**, en el orden del §2, sin tono promocional. Cierra con el código y los datos.
2. **Introduction**
   - Problema: optimización de caja negra; proliferación de metáforas sin aporte (Sörensen 2015; Aranha 2022). El NFL como contexto, no como licencia (ningún paper de CE lo usa como motivación).
   - Gap: los algoritmos bioinspirados describen comportamiento, pero ninguno traduce **leyes cuantitativas medidas** del organismo a reglas de búsqueda ni las somete a ablación (0/15 y 0/16 en el corpus).
   - Contribuciones en viñetas:
     1. tres leyes medidas de Cobetia → tres reglas estándar;
     2. ablación;
     3. evaluación con 71 funciones, 51 corridas, 3 presupuestos y estado del arte;
     4. el perfil "dónde compite y dónde no".
   - Párrafo de organización (que coincida con las secciones: falla en 4/19).
3. **Related work**, con dos ejes:
   - la taxonomía por inspiración, breve;
   - **los mecanismos equivalentes**: reducción de población (LPSR de L-SHADE), mutación de cola pesada (Cauchy y Lévy en la programación evolutiva rápida), regla 1/5 (Rechenberg) y DE de estado del arte.

   Situar cada pieza de MAREA-2 frente a su equivalente estándar (exigencia 1 de Aranha), como hacen JADE y L-SHADE con lo que heredan.
4. **De la biología a las reglas**
   - Cobetia y el paper de Frontiers.
   - **Requisitos de diseño** derivados de las tres leyes, a la manera de CMA-ES: cada ley fija una propiedad que la regla debe cumplir, y la ablación la verifica después.
   - Tabla de trazabilidad en tres columnas: **ley medida (fuente) → ecuación → componente del código**. Responde directamente a "three different things".
   - **Lista de simplificaciones** frente a la biología real, como BFO.
   - Sin vocabulario metafórico operativo.
5. **MAREA-2**
   - Algoritmo en terminología estándar, pseudocódigo, diagrama de flujo, complejidad O(·) y parámetros (`t12`).
   - Transferencia horizontal parcial: qué es en términos estándar (cruce por coordenadas con fracción adaptativa) y qué la distingue del CR adaptativo de SHADE.
   - Si hay un "padre" natural (la versión sin colonias), marcar las líneas que cambian, como jSO.
6. **Diseño experimental**
   - Suites, D, presupuestos (5k/50k/500k en FEs), 51 corridas, rivales con sus parámetros publicados, déficit de 80 FEs.
   - Métricas: media como principal, mediana como robustez, regla de 1e-8.
   - Tests: Friedman, Shaffer, Nemenyi y Wilcoxon con Holm.
   - Declarar que la ablación (15 corridas) usó las mismas 71 funciones que la evaluación: no hay separación entre diseño y evaluación (solo L-SHADE la tiene en el corpus). Una sensibilidad en otra suite lo mitigaría.
7. **Resultados**
   - Rangos (`t1`, `fig_rangos`), por benchmark (`t2`), CD (`fig_cd_nemenyi`), W/T/L (`t7`, `fig_wtl`), por clase de función (`t8`) y robustez (`t3`, `t10`).
   - Convergencia (en FEs) y boxplots, cuando los tengamos.
8. **Ablación** (`t11`): qué aporta cada ley y qué es inerte. Explicar interacciones, como JADE.
9. **Discusión**
   - Perfil de especialista: gana a 5k y en D50 5k a todos; detrás de L-SHADE/jSO a 50k/500k.
   - Familias de funciones donde pierde: híbridas y unimodales.
   - **Sección de limitaciones** explícita y con números (3/20 en el corpus).
10. **Conclusiones y trabajo futuro.**
11. **Disponibilidad de código y datos**: Zenodo.

Opcional, si hay ingeniería: además de 3–5 clásicos, una aplicación del dominio de los autores (como el diseño óptico de GWO o la celda fotovoltaica con datos medidos de KOA). En esta tesis el candidato natural es el biosensor o el control de biofilms, si cabe en el alcance.

## 9. Qué no afirmar

- "Supera al estado del arte" o "nuevo paradigma".
- "La bacteria optimiza" o "resuelve problemas".
- Términos metafóricos como nombres de operadores sin su equivalente estándar al lado.
- El resultado a 5k presentado como superioridad general.
- Un mecanismo conocido como contribución propia (el error de BKA con Cauchy).

## 10. Decisiones pendientes (del usuario y su profesor)

1. **Revista objetivo.** Swarm and Evolutionary Computation, IEEE TEVC, Evolutionary Computation, Applied Soft Computing… Decide si hacen falta problemas de ingeniería: 14 de 16 papers del género los traen y 0 de 4 de CE.
2. **Corridas nuevas en el mini:** curvas de convergencia (en FEs), tiempo, RandomSearch y sensibilidad. ¿Cuáles valen el costo?
3. **Ingeniería sí o no.** Si sí: 3–5 problemas clásicos con manejo de restricciones declarado, estadística completa, los rivales corridos por nosotros (no copiados) y la factibilidad de cada diseño verificada.
