#!/bin/bash
# Prueba cada interprete Python del sistema contra los resultados guardados
# y dice cual es el que uso el servidor Jupyter.
#
# Uso:  bash buscar_entorno.sh

echo "Buscando interpretes candidatos..."

CANDIDATOS=$(
  {
    # 1. Los que Jupyter tiene registrados como kernels
    jupyter kernelspec list --json 2>/dev/null | python3 -c "
import json, sys, os
try:
    d = json.load(sys.stdin)['kernelspecs']
except Exception:
    sys.exit()
for v in d.values():
    p = os.path.join(v['resource_dir'], 'kernel.json')
    try:
        print(json.load(open(p))['argv'][0])
    except Exception:
        pass
" 2>/dev/null
    # 2. Los del PATH
    which -a python3 python3.11 python3.12 python3.13 2>/dev/null
    # 3. Entornos conda / venv habituales
    ls -d ~/miniconda3/envs/*/bin/python3 ~/anaconda3/envs/*/bin/python3 \
          ~/miniforge3/envs/*/bin/python3 ~/.virtualenvs/*/bin/python3 \
          ~/miniconda3/bin/python3 ~/anaconda3/bin/python3 2>/dev/null
  } | sort -u
)

echo "$CANDIDATOS" | sed 's/^/  /'
echo

for PY in $CANDIDATOS; do
  [ -x "$PY" ] || continue
  echo "══════════════════════════════════════════════════════════════"
  echo "PROBANDO: $PY"
  "$PY" -c "import opfunu, numpy" 2>/dev/null || { echo "  (sin opfunu/numpy — se omite)"; echo; continue; }
  "$PY" verificar_entorno.py 2>/dev/null | grep -E "numpy|opfunu|identicas|VEREDICTO"
  # (si un interprete se cuelga, Ctrl-C pasa al siguiente)
  echo
done

echo "Busca el que diga 'ESTE ES EL ENTORNO CORRECTO'."
echo "Ese es el interprete que debes seleccionar en VS Code."
