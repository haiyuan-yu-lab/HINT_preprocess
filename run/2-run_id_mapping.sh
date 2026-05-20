SCRIPT_ROOT="/home/yl986/scripts/HINT_preprocess"

UPDATE_DIR="/home/yl986/data/HINT/update_2026"
ARCHIVE_DATA_DIR="/home/yl986/data/HINT/data/"
UPROT_SOURCE_DIR='/share/yu/resources/UniProt/release_202601/'

IDMAPPING_DIR=${UPROT_SOURCE_DIR}/"knowledgebase/idmapping/"
UNMAPPED_LOG="unmapped_2026.log"
cd $SCRIPT_ROOT
python -u 3-create_idmapping.py --update_dir ${UPDATE_DIR} --id_mapping_source_dir ${IDMAPPING_DIR} --unmapped_log_file ${UNMAPPED_LOG}
