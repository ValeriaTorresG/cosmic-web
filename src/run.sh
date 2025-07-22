#!/usr/bin/env bash
set -e

# Sitúate en la raíz del repo
cd "$(dirname "$0")"/..

# Rutas fijas
REGION_DIR="/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/regions"
ENTROPY_DIR="/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/entropy"
CATALOG_DIR="/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/catalog"

# Asegura que existan
mkdir -p "$ENTROPY_DIR" "$CATALOG_DIR"

echo "1) Generando posiciones, conteos, pares y clasificación (entropía)…"
# Este script ya sabe leer de REGION_DIR y escribir en ENTROPY_DIR internamente
python src/entropy.py

echo
echo "2) Construyendo catálogos friends-of-friends…"
# web_catalog.py espera solo base_dir, out_dir, webtype, void_limit, knot_limit
python src/web_catalog.py \
    --base_dir   "$ENTROPY_DIR" \
    --out_dir    "$CATALOG_DIR" \
    --webtype    all \
    --void_limit -0.9 \
    --knot_limit 0.9

echo
echo "✅ Pipeline completado"
echo "  • Subregiones filtradas: $REGION_DIR"
echo "  • Resultados entropía:   $ENTROPY_DIR"
echo "  • Catálogos finales:     $CATALOG_DIR"
