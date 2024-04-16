UNIPARC_DIR='/home/yl986/data/HINT/uniprot_source/release_202401/uniparc/'
OUTPUT_DIR='/home/yl986/data/HINT/outputs_2023/'

for fname in $(ls $UNIPARC_DIR | grep "xml.gz")
do
    python parse_uniparc.py $UNIPARC_DIR $fname $OUTPUT_DIR
done
