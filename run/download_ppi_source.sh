#!/usr/bin/bash
# set -euo pipefail

# UPDATE_DIR="/home/yl986/data/HINT/update_2026"
# ARCHIVE_DATA_DIR="/home/yl986/data/HINT/data/"
UPDATE_DIR=$1
ARCHIVE_DATA_DIR=$2

DATA_ROOT=${UPDATE_DIR}/"data/"

# WIPE="${WIPE:-0}"  # set WIPE=1 to rm -rf UPDATE_DIR
# [ ${WIPE} -eq 1 ] && rm -rf ${UPDATE_DIR}
if [ ! -d ${DATA_ROOT} ]; then
    mkdir -p ${DATA_ROOT}
fi
echo "Prepare PPI data source in ${DATA_ROOT}"

# Copy required inputs
echo "# Copy static datasets..."
static_fnames=$(ls ${ARCHIVE_DATA_DIR}/"static_datasets")
cp -a ${ARCHIVE_DATA_DIR}/"static_datasets" ${DATA_ROOT}
for file in $static_fnames; do
    echo " - $(basename ${file}) copied"
done
echo "# Copy bouncer (manually curated interaction whitelist / blacklist)..."
cp -a ${ARCHIVE_DATA_DIR}/"bouncer" ${DATA_ROOT}
echo "# Copy explicit_interactomes (manually reviewed)..."
cp -a ${ARCHIVE_DATA_DIR}/"explicit_interactomes" ${DATA_ROOT}
echo "# Copy taxid2name.txt..."
cp -p ${ARCHIVE_DATA_DIR}/"taxid2name.txt" ${DATA_ROOT}
echo "# Copy evidence code reference..."
cp -p ${ARCHIVE_DATA_DIR}/evidence_code*.txt ${DATA_ROOT}

# Create required output structure
mkdir -p ${UPDATE_DIR}/"outputs/HINT_format/taxa"
mkdir -p ${UPDATE_DIR}/"caches"

# Downloads
echo "# Download active data sources (BioGrid, IntAct)"
BIOGRID_URL="https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.mitab.zip"
INTACT_URL="ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip"

TMP=${DATA_ROOT}/tmp
PARSE=${DATA_ROOT}/parseTargets
STATIC=${DATA_ROOT}/static_datasets

mkdir -p ${STATIC}

# Clean and recreate working dirs
rm -rf ${PARSE} ${TMP}
mkdir -p ${PARSE} ${TMP}

# BioGRID
wget -P ${PARSE} -q ${BIOGRID_URL}
unzip -d ${PARSE} ${PARSE}/${BIOGRID_URL##*/}

# IntAct
wget -P ${PARSE} -q ${INTACT_URL}
unzip -d ${PARSE} ${PARSE}/${INTACT_URL##*/}

# Copy static datasets
cp -p ${STATIC}/* ${PARSE}/

# Copy IRES files
IRES_FPATH="/share/yu/resources/local_resources/pdb/parsed_file/ires_perpdb_w_sifts_mapping.txt"
echo "# Copy curated PDB information (IRES, pdb_info.txt)"
cp -p ${IRES_FPATH} ${PARSE}/
PDB_INFO_PATH="/home/yl986/data/HINT/pdb_parsed"
cp -p ${PDB_INFO_PATH}/* ${PARSE}

# Cleanup
rm -rf ${TMP}

echo "Finished downloading PPI source data."