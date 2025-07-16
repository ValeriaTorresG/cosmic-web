set -e

cd "$(dirname "$0")"/..

echo "Running pos,count,pairs"
python src/entropy.py

RESULTS_DIR="data/results"
CATALOGS_DIR="data/catalogs"
mkdir -p "${CATALOGS_DIR}"

echo "----- Generating catalogs"
for rosette in $(seq 0 19); do
  echo "-- Rosette ${rosette}"
  python src/make_web_catalog.py \
    --posfile     "${RESULTS_DIR}/rosette_${rosette}_pos.txt" \
    --pairfile    "${RESULTS_DIR}/rosette_${rosette}_pairs.txt" \
    --countfile   "${RESULTS_DIR}/rosette_${rosette}_counts.txt" \
    --webtype      void \
    --catalogfile "${CATALOGS_DIR}/rosette_${rosette}_catalog_all.csv"
done

echo "----- Funcionaaaa, guarda en ${CATALOGS_DIR}"
