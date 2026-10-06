# Ficha de extracción (una por paper)

Lee el PDF COMPLETO (texto con `python3 -c "from pypdf import PdfReader ..."`; si hace falta, mira figuras/tablas con la herramienta Read sobre el PDF y `pages`).
Escribe `<scratchpad>/corpus/NN.json` (NN = número del archivo) con EXACTAMENTE estas claves. Valores: true/false/null (null = no aplica o no se puede determinar) más un campo `_nota` corto cuando haga falta matizar. Todo texto libre en español, breve. Cita página cuando afirmes un defecto.

{
 "id": "NN", "algoritmo": "", "autores": "", "revista": "", "anio": 0, "doi": "", "paginas": 0,
 "tipo": "metafora | sin_metafora | microbiano",
 "inspiracion": "",                       // fuente de la analogía, o "ninguna"
 "secciones": [],                         // títulos de sección de primer nivel, en orden
 "abstract_orden": "",                    // secuencia de movimientos del abstract, p.ej. "propuesta→fases→validación→resultado→ingeniería→código"
 "nfl_como_motivacion": bool,
 "taxonomia_familias": bool, "tabla_ventajas_desventajas": bool,
 "revision_literatura_seccion_propia": bool,
 "contribuciones_en_lista": bool,
 "parrafo_organizacion": bool, "organizacion_coincide": bool,
 "seccion_inspiracion_con_figuras": bool,
 "metafora_modelada_en_ecuaciones": "total | parcial | no | n/a",   // ¿cada elemento de la historia tiene ecuación y viceversa?
 "datos_biologicos_medidos_en_ecuaciones": bool,
 "operadores_prestados": "",              // mecanismos conocidos (mutación DE, Lévy, greedy, LPSR...) y si los DECLARA como tales
 "fases_exploracion_explotacion_explicitas": bool, "exploracion_medida": bool,
 "pseudocodigo": bool, "diagrama_flujo": bool,
 "complejidad_O": bool, "complejidad_correcta": "sí | no | dudosa | n/a",  // justifica si no
 "suites": [],                            // p.ej. ["clásicas F1-F23", "CEC2017"], con dimensiones
 "dimensiones": [],
 "presupuesto": "",                       // máx FEs o iteraciones×población tal como lo dan
 "compara_por": "FEs | iteraciones | ambiguo",
 "fes_por_agente_iter_del_propuesto": "",// cuántas evaluaciones gasta por agente e iteración, si se deduce
 "presupuesto_nivel_cec_10000D": bool,
 "corridas": 0,
 "competidores": [], "n_competidores": 0,
 "incluye_LSHADE": bool, "incluye_jSO": bool, "incluye_CMAES": bool, "incluye_RandomSearch": bool,
 "parametros_rivales_publicados": bool,
 "tests": [],                             // Wilcoxon rank-sum/signed-rank, Friedman, ...
 "posthoc": "",                           // Holm, Nemenyi, Bonferroni-Dunn, ninguno
 "correccion_multiples": bool,
 "conteo_WTL": bool,
 "curvas_convergencia": bool, "convergencia_eje": "FEs | iteraciones | n/a",
 "analisis_cualitativo_4_paneles": bool,
 "boxplots": bool,
 "diversidad_porcentaje_expl": bool,
 "sensibilidad_parametros": bool, "sensibilidad_en_mismas_funciones": bool,
 "escalabilidad_dimension": bool,
 "tiempo_computo": bool,
 "ablacion_por_componente": bool,
 "problemas_ingenieria": [],              // nombres
 "problemas_ingenieria_reales": bool,     // aplicación real además de los de manual
 "ingenieria_rivales_copiados_literatura": bool,
 "ingenieria_infactibles_detectadas": "", // comprueba el mejor diseño reportado de 1–2 problemas clásicos (p.ej. recipiente a presión: x1,x2 múltiplos de 0.0625; costo no menor que ~5885.33 continuo / 6059.71 discreto; viga soldada ~1.6702/1.7249 según versión; resorte ~0.012665) y di qué encontraste
 "codigo_publico": bool, "datos_crudos_publicos": bool,
 "seccion_limitaciones": "honesta | genérica | no",
 "reporta_donde_pierde": bool,
 "inconsistencias": "",                   // nº de rivales/problemas distinto entre abstract y texto, valores imposibles, erratas importantes
 "estilo_escritura": "",                  // 2–3 frases: tono, longitud, uso de vocabulario metafórico, cómo presenta resultados
 "lo_mejor_para_copiar": "",              // 1–2 cosas concretas que un paper riguroso debería imitar de este
 "_nota": ""
}

Para los papers sin metáfora (JADE, L-SHADE, jSO, CMA-ES) los campos de metáfora van "n/a"/null; fíjate especialmente en CÓMO justifican el diseño sin metáfora (análisis, derivación, ablación, comparación con el padre).
