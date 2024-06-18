# Run IRES calculation pipeline

# Add directory of libgfortran.so.4 to LD_LIBRARY_PATH to run `naccess` program
# export LD_LIBRARY_PATH='/home/yl986/bin/miniconda3/envs/masif/lib/':$LD_LIBRARY_PATH

IRES_PARSED_DIR='/local/storage/yl986/script_hub/HINT_preprocess/pdb_data_prep/data/ires/parsed_files'
dist_file_prefix='bounding_box_pdb_chain_dists'

python ./pdb_data_prep/split_to_batch.py \
    -input_dir $IRES_PARSED_DIR \
    -output_dir $IRES_PARSED_DIR/batch_inputs \
    -fname $dist_file_prefix'.txt' \
    -chunksize 100000

n_batch=$(ls $IRES_PARSED_DIR/batch_inputs | grep -c $dist_file_prefix)

for (( IDX=0; IDX<n_batch; IDX++ ))
do
# echo $IDX
    nohup python -u ./pdb_data_prep/2-ires_cal_yl.py \
        -fdist batch_inputs/bounding_box_pdb_chain_dists_0_$IDX'.txt' \
        -output ires_perpdb_alltax_0_$IDX'.txt' \
        -fneg NOires_perpdb_alltax.txt >logs/out0_$IDX 2>logs/err0_$IDX &
done
        # 

# nohup python -u ./pdb_data_prep/ires_for_pdb_bundle.py >logs/pdblike.o 2>logs/pdblike.e