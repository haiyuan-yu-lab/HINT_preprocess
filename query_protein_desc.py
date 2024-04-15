import os, sys
from pathlib import Path

import pandas as pd


"""
Extract description for selected protein IDs from meta info file (extracted from FASTA)
    Input:  protein-level ID mapping (source ID -> UniProt ID) (`id_to_uniprot.txt`)
    Output: description table for selected proteins (`uniprot_descriptions.txt`)
"""

def extract_uprot_meta_info(df_prot, source_path, out_path, pid_col='UniProt', meta_pid_col='UniProt', 
                            meta_prefix=['sprot', 'sprot_varsplic', 'trembl'], chunksize=1e6, output_name='uniprot_descriptions.txt'):
    """
    Extract description for selected protein IDs
    """
    if isinstance(source_path, str):
        source_path = Path(source_path)

    if isinstance(out_path, str):
        out_path = Path(out_path)
    
    if not out_path.exists():
        out_path.mkdir(parents=True)

    prot_remain = set(df_prot[pid_col])
    print('# UniProt IDs:', len(prot_remain))

    header = True
    mode = 'w'
    for prefix in meta_prefix:
        fname = f'{prefix}_meta.txt'
        if not (source_path / fname).exists():
            print('{} not found!'.format(str(source_path / fname)))
            continue

        chunks = pd.read_csv(source_path / fname, sep='\t', chunksize=chunksize)
        for df_meta in chunks:
            records = df_meta[df_meta[meta_pid_col].isin(prot_remain)]
            records.to_csv(out_path / output_name, mode=mode, header=header, sep='\t', index=False)
            prot_remain = prot_remain - set(records[meta_pid_col])
            print('# UniProt IDs remain:', len(prot_remain))
            if len(prot_remain) == 0:
                break
            header = False  # write header only the first time
            mode = 'a'


if __name__ == '__main__':
    data_root = Path('/home/yl986/data/HINT/outputs_2023/cache')
    prot_meta_root = Path('/home/yl986/data/HINT/uniprot_source/release_202401/knowledgebase/complete/meta')
    df_prot = pd.read_csv(data_root / 'id_to_uniprot.txt', sep='\t')

    extract_uprot_meta_info(df_prot, prot_meta_root, out_path=data_root, pid_col='primary_ac_short')