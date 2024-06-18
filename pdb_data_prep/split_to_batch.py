import os, sys
import argparse
import pandas as pd

if __name__ == '__main__':
    ## I/O arguments at command line
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Split target file into small batches')
    parser.add_argument('-input_dir', help='directory of input file')
    parser.add_argument('-output_dir', default='', help='directory of output file (small batches)')
    parser.add_argument('-fname', default='bounding_box_pdb_chain_dists.txt', help='target file')
    parser.add_argument('-chunksize', default=1e4, type=int, help='Chunksize of each batch file')

    args = parser.parse_args()

    data_dir = args.input_dir
    output_dir = args.output_dir
    if not output_dir:
        output_dir = os.path.join(data_dir, 'batch_inputs')

    fname = args.fname
    chunksize = args.chunksize

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    prefix = fname.split('.')[0]

    chunks = pd.read_csv(os.path.join(output_dir, fname), sep='\t', chunksize=chunksize)
    
    for i, df_chunk in enumerate(chunks):
        df_chunk.to_csv(os.path.join(output_dir, f'{prefix}_{i}.txt'), index=False, sep='\t')
