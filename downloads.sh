# Download data from UniProt FTP site (current release)

# SOURCE_DIR='/home/yl986/data/HINT/uniprot_source/release_202401/'
# SOURCE_DIR='/share/yu/resources/UniProt/release_202601'
SOURCE_DIR=$1
UNIPROT_FTP='https://ftp.uniprot.org/pub/databases/uniprot/current_release/'

TMP_FLIST="$SOURCE_DIR/tmp_flist"

if [ ! -d "$SOURCE_DIR" ]; then
    mkdir -p "$SOURCE_DIR"
fi

list_remote_files () {
    local url="$1"
    local suffix="$2"
    curl -fsSL "$url" \
      | sed -n 's/.*href="\([^"]*\)".*/\1/p' \
      | grep -E "^[-A-Za-z0-9._]+${suffix}$" \
      | sort -u
}

# Download id mapping files
base_url=$UNIPROT_FTP'knowledgebase/idmapping/'
suffix='dat.gz'

mkdir -p $SOURCE_DIR/knowledgebase/idmapping
echo "Downloading id mapping files..."

list_remote_files "$base_url" "$suffix" > "$TMP_FLIST"
while IFS= read -r i; do
    wget -nc "${base_url}${i}" -P "$SOURCE_DIR/knowledgebase/idmapping/"
done < "$TMP_FLIST"


# Download uniprot docs
base_url=$UNIPROT_FTP'knowledgebase/complete/docs/'
suffix='txt'

mkdir -p $SOURCE_DIR/knowledgebase/docs
echo "Downloading UniProt documentation files..."

list_remote_files "$base_url" "$suffix" > "$TMP_FLIST"
while IFS= read -r i; do
    wget -nc "${base_url}${i}" -P "$SOURCE_DIR/knowledgebase/docs/"
done < "$TMP_FLIST"

# Downlaod uniprot fasta files
# uniprot_sprot.fasta.gz, uniprot_sprot_varsplic.fasta.gz, uniprot_trembl.fasta.gz
base_url=$UNIPROT_FTP'knowledgebase/complete/'
suffix='fasta.gz'

mkdir -p $SOURCE_DIR/knowledgebase/complete
echo "Downloading UniProt FASTA files..."

list_remote_files "$base_url" "$suffix" > "$TMP_FLIST"
while IFS= read -r i; do
    wget -nc "${base_url}${i}" -P "$SOURCE_DIR/knowledgebase/complete/"
done < "$TMP_FLIST"


# Download UniParc XML files
base_url=$UNIPROT_FTP'uniparc/xml/all/'
#wget --no-directories --no-verbose --no-remove-listing -O - $base_url
suffix='xml.gz'

mkdir -p $SOURCE_DIR/uniparc
echo "Downloading UniParc XML files..."

list_remote_files "$base_url" "$suffix" > "$TMP_FLIST"
while IFS= read -r i; do
    wget -nc "${base_url}${i}" -P "$SOURCE_DIR/uniparc"
done < "$TMP_FLIST"

echo "UniProt knowledgebase download complete!"