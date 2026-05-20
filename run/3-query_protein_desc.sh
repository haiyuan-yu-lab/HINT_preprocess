SCRIPT_ROOT="/home/yl986/scripts/HINT_preprocess"

UPDATE_DIR="/home/yl986/data/HINT/update_2026"
ARCHIVE_DATA_DIR="/home/yl986/data/HINT/data/"
UPROT_SOURCE_DIR='/share/yu/resources/UniProt/release_202601/'

UNMAPPED_LOG="unmapped_2026.log"
cd $SCRIPT_ROOT
# Pre-requisite: `id_to_uniprot.txt` file created (`4-HINT_data_curation.ipynb`)
python -u query_protein_desc.py --update_dir ${UPDATE_DIR} --uniprot_knowledgebase ${UPROT_SOURCE_DIR}/"knowledgebase/"
