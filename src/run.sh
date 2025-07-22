#!/usr/bin/env bash
set -e

# root of folder
cd "$(dirname "$0")"/..

REGION_DIR="/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/regions"
ENTROPY_DIR="/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/entropy"
CATALOG_DIR="/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/catalog"

mkdir -p "$ENTROPY_DIR" "$CATALOG_DIR"

echo "Running pre processing, entropy script"
python src/entropy.py

echo
echo "Constructing the web catalog"
python src/web_catalog.py \
    --base_dir   "$ENTROPY_DIR" \
    --out_dir    "$CATALOG_DIR" \
    --webtype    all \
    --void_limit -0.9 \
    --knot_limit 0.9

echo