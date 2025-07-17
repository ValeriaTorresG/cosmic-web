set -e

cd "$(dirname "$0")"/..

echo "Running pos,count,pairs"
<<<<<<< HEAD
python src/entropy.py
=======
python src/entropy_c.py
>>>>>>> origin/main

RESULTS_DIR="data/results"
CATALOGS_DIR="data/catalogs"
mkdir -p "${CATALOGS_DIR}"

echo "----- Generating catalogs"
for rosette in $(seq 0 19); do
  echo "-- Rosette ${rosette}"
  python src/catalog.py \
    --posfile     "${RESULTS_DIR}/rosette_${rosette}_pos.txt" \
    --pairfile    "${RESULTS_DIR}/rosette_${rosette}_pairs.txt" \
    --countfile   "${RESULTS_DIR}/rosette_${rosette}_counts.txt" \
<<<<<<< HEAD
    --webtype      filament \
    --catalogfile "${CATALOGS_DIR}/rosette_filament_${rosette}_catalog_all.csv"
=======
    --webtype      void \
    --catalogfile "${CATALOGS_DIR}/rosette_${rosette}_catalog_all.csv"
>>>>>>> origin/main
done

echo "----- Funcionaaaa, guarda en ${CATALOGS_DIR}"
