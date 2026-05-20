#!/bin/bash
#SBATCH -J 0_prep_uprot_source_                         # Job name
#SBATCH -o logs/0_prep_uprot_source_%j.out                  # output file (%j expands to jobID)
#SBATCH -e logs/0_prep_uprot_source_%j.err                  # error log file (%j expands to jobID)
#SBATCH --mail-user=yl986@cornell.edu        # Email address to send results to.
#SBATCH -N 1                                 # Total number of nodes requested
#SBATCH -n 4                                 # Total number of cores requested
#SBATCH --get-user-env                       # retrieve the users login environment
#SBATCH --mem=15G                           # server memory requested (per node)
#SBATCH -t 2-00:00:00                           # Time limit (hh:mm:ss)
#SBATCH --partition=yu                       # Request partition

TARGET_DIR='/share/yu/resources/UniProt/release_202601'
DOWNLOAD=FALSE
PARSE=TRUE

if [ "$DOWNLOAD" = TRUE ]; then
    sh download_uniprot.sh $TARGET_DIR
fi

if [ "$PARSE" = TRUE ]; then
    source activate py38
    echo "Parsing UniProt source data..."
    python -u ../1-parse_source_data.py ${TARGET_DIR}/"knowledgebase"
fi