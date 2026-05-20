#!/bin/bash
#SBATCH -J hint_desc_query                         # Job name
#SBATCH -o hint_desc_query_%j.out                  # output file (%j expands to jobID)
#SBATCH -e hint_desc_query_%j.err                  # error log file (%j expands to jobID)
#SBATCH --mail-user=yl986@cornell.edu        # Email address to send results to.
#SBATCH -N 1                                 # Total number of nodes requested
#SBATCH -n 4                                 # Total number of cores requested
#SBATCH --get-user-env                       # retrieve the users login environment
#SBATCH --mem=40G                           # server memory requested (per node)
#SBATCH -t 2-00:00:00                           # Time limit (hh:mm:ss)
#SBATCH --partition=yu                       # Request partition


# python -u query_protein_desc.py
source activate py38
sh 3-query_protein_desc.sh
