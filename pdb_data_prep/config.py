OUTPUT_DIR = '/local/storage/yl986/script_hub/HINT_preprocess/pdb_data_prep/data/'
# PDB_DATA_DIR = '/local/storage/resources/pdb/data/'
# PDB_DATA_DIR = '/local/storage/yl986_02/data/pdb/'  # server02
PDB_DATA_DIR = '/fs/cbsuhy02/storage/yl986_02/data/pdb/'  # server01
PDB_BUNDLE_DATA_DIR = f'{PDB_DATA_DIR}/pdb_like/data/'
# # PDB_TO_RUN = '/local/storage/yl986/script_hub/nightly_scripts/test_data/pdb_to_run.txt'
PDB_TO_RUN = '/local/storage/yl986/script_hub/HINT_preprocess/pdb_data_prep/data/pdb_to_run.txt'
PDB_BUNDLE_LIST = '/local/storage/yl986/script_hub/HINT_preprocess/pdb_data_prep/data/pdb/pdb_like/parsed_files/pdb_like_list.txt'
# g2 configs
# OUTPUT_DIR = '/home/yl986/scripts/HINT_preprocess/pdb_data_prep/data/test_data'
#PDB_DATA_DIR = '/share/yu/resources/pdb/'
# PDB_TO_RUN = '/home/yl986/scripts/HINT_preprocess/pdb_data_prep/data/pdb_to_run.txt'

# SIFTS: https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
# SIFTS_DATA_DIR = '/local/storage/resources/sifts/data/'
# SIFTS_DATA_ROOT = '/local/storage/yl986/script_hub/nightly_scripts/test_data/sifts2024/'
SIFTS_DATA_ROOT = '/local/storage/resources/sifts/20240603'
SIFTS_DATA_DIR = f'{SIFTS_DATA_ROOT}/data/'
# SIFTS_PARSED_DIR = f'{SIFTS_DATA_ROOT}/parsed_files/'
SIFTS_PARSED_DIR = f'{OUTPUT_DIR}/sifts/parsed_files'
SIFTS_MAPPING_FILE = f'{SIFTS_PARSED_DIR}/pdbresiduemapping.txt'


# IRES: interface residues
IRES_PARSED_DIR = f'{OUTPUT_DIR}/ires/parsed_files/'  # To-Do: modify
IRES_NEG_DIR = f'{IRES_PARSED_DIR}/_negative_results'

DIST_THRES = 8
NACCESS_PATH = '/local/storage/yl986/script_hub/naccess2.1.1/naccess'