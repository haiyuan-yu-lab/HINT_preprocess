#!/bin/bash
#SBATCH -J pdb_prep_                         # Job name
#SBATCH -o logs/pdb_prep_%j.out                  # output file (%j expands to jobID)
#SBATCH -e logs/pdb_prep_%j.err                  # error log file (%j expands to jobID)
#SBATCH --mail-user=yl986@cornell.edu        # Email address to send results to.
#SBATCH -N 1                                 # Total number of nodes requested
#SBATCH -n 4                                 # Total number of cores requested
#SBATCH --get-user-env                       # retrieve the users login environment
#SBATCH --mem=15G                           # server memory requested (per node)
#SBATCH -t 2-00:00:00                           # Time limit (hh:mm:ss)
#SBATCH --partition=yu                       # Request partition

OUTPUT_DIR="/home/yl986/data/HINT/pdb_parsed"

source activate py38
python -u ../00-extract_pdb_info.py --output_root $OUTPUT_DIR