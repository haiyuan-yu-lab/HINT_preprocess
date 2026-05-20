#!/bin/bash
#SBATCH -J 2_id_mapping_                         # Job name
#SBATCH -o logs/2_id_mapping_%j.out                  # output file (%j expands to jobID)
#SBATCH -e logs/2_id_mapping_%j.err                  # error log file (%j expands to jobID)
#SBATCH --mail-user=yl986@cornell.edu        # Email address to send results to.
#SBATCH -N 1                                 # Total number of nodes requested
#SBATCH -n 4                                 # Total number of cores requested
#SBATCH --get-user-env                       # retrieve the users login environment
#SBATCH --mem=50G                           # server memory requested (per node)
#SBATCH -t 2-00:00:00                           # Time limit (hh:mm:ss)
#SBATCH --partition=yu                       # Request partition

source activate py38

sh 2-run_id_mapping.sh
