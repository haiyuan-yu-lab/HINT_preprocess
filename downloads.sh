# Download data from UniProt FTP site (current release)

SOURCE_DIR='/home/yl986/data/HINT/uniprot_source/release_202401/'
UNIPROT_FTP='https://ftp.uniprot.org/pub/databases/uniprot/current_release/'

# Download id mapping files
base_url=$UNIPROT_FTP'knowledgebase/idmapping/'
suffix='dat.gz'

mkdir -p $SOURCE_DIR/knowledgebase/idmapping

curl --silent $base_url | grep -o 'href=".*">' | sed 's/href="//;s/\/">//' | sed 's/">//g' | grep $suffix > ./tmp_flist

for i in $(cat ./tmp_flist)
do
    wget -nc $base_url/$i -P $SOURCE_DIR/knowledgebase/idmapping/
done


# Download uniprot docs
base_url=$UNIPROT_FTP'knowledgebase/complete/docs/'
suffix='txt'

mkdir -p $SOURCE_DIR/knowledgebase/docs

curl --silent $base_url | grep -o 'href=".*">' | sed 's/href="//;s/\/">//' | sed 's/">//g' | grep $suffix > ./tmp_flist

for i in $(cat ./tmp_flist)
do
    wget -nc $base_url/$i -P $SOURCE_DIR/knowledgebase/docs/
done


# Downlaod uniprot fasta files
# uniprot_sprot.fasta.gz, uniprot_sprot_varsplic.fasta.gz, uniprot_trembl.fasta.gz
base_url=$UNIPROT_FTP'knowledgebase/complete/'
suffix='fasta.gz'

mkdir -p $SOURCE_DIR/knowledgebase/complete

curl --silent $base_url | grep -o 'href=".*">' | sed 's/href="//;s/\/">//' | sed 's/">//g' | grep $suffix > ./tmp_flist

for i in $(cat ./tmp_flist)
do
    wget -nc $base_url/$i -P $SOURCE_DIR/knowledgebase/complete/
done


# Download UniParc XML files
base_url=$UNIPROT_FTP'uniparc/xml/all/'
#wget --no-directories --no-verbose --no-remove-listing -O - $base_url
suffix='xml.gz'

mkdir -p $SOURCE_DIR/uniparc

curl --silent $base_url | grep -o 'href=".*">' | sed 's/href="//;s/\/">//' | sed 's/">//g' | grep $suffix > ./tmp_flist

for i in $(cat ./tmp_flist)
do
    wget -nc $base_url/$i -P $SOURCE_DIR/uniparc
done


rm ./tmp_flist