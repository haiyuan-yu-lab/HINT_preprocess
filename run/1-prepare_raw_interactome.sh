SCRIPT_ROOT="/home/yl986/scripts/HINT_preprocess"

UPDATE_DIR="/home/yl986/data/HINT/update_2026"
ARCHIVE_DATA_DIR="/home/yl986/data/HINT/data/"
DOWNLOAD=FALSE

if [ $DOWNLOAD = TRUE ]; then
  sh download_ppi_source.sh $UPDATE_DIR $ARCHIVE_DATA_DIR
fi

# Prerequisite: download PPI source data (`download_ppi_source.sh`)
cd $SCRIPT_ROOT
python -u 2-prepare_dataset.py --update_dir $UPDATE_DIR --archive_data_dir $ARCHIVE_DATA_DIR