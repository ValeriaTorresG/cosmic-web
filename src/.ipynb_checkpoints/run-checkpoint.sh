set -e

cd "$(dirname "$0")"/..

echo "Running pos,count,pairs"
python src/entropy.py


RESULTS_DIR="/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/results"
CATALOGS_DIR="/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/catalogs"
mkdir -p "${CATALOGS_DIR}"

# echo "----- Generating catalogs"
# for rosette in $(seq 0 19); do
#   echo "-- Rosette ${rosette}"
#   python src/catalog.py \
#     --posfile     "${RESULTS_DIR}/rosette_${rosette}_pos.txt" \
#     --pairfile    "${RESULTS_DIR}/rosette_${rosette}_pairs.txt" \
#     --countfile   "${RESULTS_DIR}/rosette_${rosette}_counts.txt" \
#     --webtype      filament \
#     --catalogfile "${CATALOGS_DIR}/rosette_filament_${rosette}_catalog_all.csv"
# done

echo "----- Funcionaaaa, guarda en ${CATALOGS_DIR}"
