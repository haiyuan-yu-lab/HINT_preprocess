import os, sys
import json
import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path
import time
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from constants import ACCEPTED_IDS
# from configs import *

def fetch_select_type_to_uprot(fpath, fout_path, target_ids, chunksize=1000000):
    """
    Generate ID mapping reference from `idmapping.dat.gz` (from UniProt FTP site) 
    and write source-to-uniprot ID mapping to JSON file

    Args:
        fpath: path to idmapping.dat.gz file
        fout_path: output file path
        target_ids: target ID keywords
        chunksize: number of records in each chunk [default: 1000000]
    """
    idmap_by_type = defaultdict(dict)
    chunks = pd.read_csv(fpath, sep='\t', names=['uprot', 'id_type', 'id'], chunksize=chunksize)

    start_time = time.time()
    n_processed = 0
    for df in chunks:
        n_processed += len(df)
        print(f"Processed {n_processed} records... (Time elapsed: {time.time() - start_time:.2f}s)")
        df_tmp = df[df['id_type'].isin(target_ids)]
        for itype in target_ids:
            idx = (df_tmp['id_type'] == itype)
            agg_id = df_tmp.loc[idx].groupby('id').agg({'uprot': lambda x: '|'.join(sorted(set(x)))}).reset_index()
            
            cur_dict = dict(zip(agg_id['id'], agg_id['uprot']))
            for cur_id, cur_uprot in cur_dict.items():
                if cur_id in idmap_by_type[itype]:
                    uprot_lst = idmap_by_type[itype][cur_id].split('|') + cur_uprot.split('|')
                    idmap_by_type[itype][cur_id] = '|'.join(set(uprot_lst)).strip('|')
                else:
                    idmap_by_type[itype][cur_id] = cur_uprot
            
    print(f"\nFinished processing {n_processed} records in {time.time() - start_time:.2f}s")
    with open(fout_path, 'w') as f:
        json.dump(idmap_by_type, f, indent=2)
    
    return idmap_by_type


def fetch_protein_gene_desc(fpath, fout_path, target_fields, chunksize=1000000):
    """
    Fetch protein description from `idmapping.dat.gz` (from UniProt FTP site) 
    and write protein-to-description mapping to JSON file
    Args:
        fpath: path to idmapping.dat.gz file
        fout_path: output file path
        target_fields: target information keywords to look for
        chunksize: number of records in each chunk [default: 1000000]
    
    Returns:
        dictionary maps UniProt ID to information dict
    """

    chunks = pd.read_csv(fpath, sep='\t', names=['uprot', 'id_type', 'id'], chunksize=chunksize)
    mode = 'w'
    output_header = True
    output_cols = ['uprot'] + target_fields
    uprot2info_dict = defaultdict(dict)

    start_time = time.time()
    n_processed = 0
    for df in chunks:
        n_processed += len(df)
        print(f"Processed {n_processed} records... (Time elapsed: {time.time() - start_time:.2f}s)")
        df_info = df[df['id_type'].isin(target_fields)].astype({'id': str}).\
            groupby(['uprot', 'id_type']).agg({'id': lambda x: '|'.join(sorted(set(x)))}).reset_index()
        df_info = df_info.pivot(columns=['id_type'], values='id', index='uprot').fillna('')
        for field in target_fields:
            if field not in df_info.columns:
                df_info[field] = ''
        cur_dict = df_info.to_dict(orient='index')
        for prot, cur_desc in cur_dict.items():
            if prot in uprot2info_dict:
                for field in target_fields:
                    if cur_desc[field] != '':
                        if field in ['Gene_Name', 'Gene_ORFName']:
                            id_lst = uprot2info_dict[prot][field].split('|') + cur_desc[field].split('|')
                            uprot2info_dict[prot][field] = '|'.join(id_lst).strip('|')
                        else:
                            uprot2info_dict[prot][field] = cur_desc[field]
            else:
                uprot2info_dict[prot] = cur_desc
    
    print(f"\nFinished processing {n_processed} records in {time.time() - start_time:.2f}s")
    with open(fout_path, 'w') as f:
        json.dump(uprot2info_dict, f, indent=2)

    return uprot2info_dict

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Create ID mapping')
    parser.add_argument('--update_dir', type=str, help='Path to the update directory')
    parser.add_argument('--id_mapping_source_dir', type=str, default='/share/yu/resources/UniProt/release_202601/knowledgebase/idmapping/', 
                        help='Directory for UniProt id-mapping file')
    parser.add_argument('--unmapped_log_file', type=str, default='unmapped.log', help='Records of nmapped protein IDs')
    parser.add_argument('--overwrite_source_cache', action='store_true', help='Whether to overwrite existing mapping files')
    args = parser.parse_args()
    source_root = Path(args.id_mapping_source_dir)
    update_root = Path(args.update_dir)
    # source_root = Path(ID_MAPPING_SOURCE_DIR)
    # update_root = Path('/home/yl986/data/HINT')
    # update_root = Path(UPDATE_DIR)

    idmapping_source_fpath = source_root / 'idmapping.dat.gz'
    if not idmapping_source_fpath.exists():
        print(f'ID mapping source file not found: {idmapping_source_fpath}')
        exit(1)

    idmapping_out_fpath = source_root / 'cache/target_type_to_uprot.json'
    prot_to_desc_json = source_root / 'cache/prot_gene_info.json'
    prot_to_desc_tab = source_root / 'cache/prot_gene_info.tsv'
    if not (source_root / 'cache').exists():
        (source_root / 'cache').mkdir(parents=True)

    # OVERWRITE_SOURCE_MAPPING = False
    OVERWRITE_SOURCE_MAPPING = args.overwrite_source_cache

    target_ids = ['BioGRID',
                'ChEMBL',
                'ComplexPortal',
                'DIP',
                'EMBL',
                'Ensembl',
                'EnsemblGenome',
                'GeneID',
                'PDB',
                'Reactome',
                'RefSeq']
    
    print('Fetching protein ID mapping...')

    if not idmapping_out_fpath.exists() or OVERWRITE_SOURCE_MAPPING:
        print(f'Parsing all protein ID mapping from {idmapping_source_fpath}...')
        idmap_by_type = fetch_select_type_to_uprot(idmapping_source_fpath, 
                                                   idmapping_out_fpath,
                                                   target_ids=target_ids,
                                                   chunksize=1000000)
        print(f'Saved protein ID mapping to {idmapping_out_fpath}')
    else:
        with open(idmapping_out_fpath, 'r') as f:
            idmap_by_type = json.load(f)

    with open(update_root / 'outputs/cache/mapping_targets_by_type.json', 'r') as f:
        mapping_targets_by_type = json.load(f)

    result_dict = dict()
    unmapped_ids = []
    total_counts = 0
    for itype, target_lst in mapping_targets_by_type.items():
        if itype not in ACCEPTED_IDS:
            continue
        itype_key = ACCEPTED_IDS[itype]
        total_counts += len(target_lst)
        for target_id in target_lst:
            try:
                result_dict[f'{itype_key}|{target_id}'] = idmap_by_type[itype_key][target_id]
            except KeyError:
                try:
                    target_id = target_id.upper()
                    result_dict[f'{itype_key}|{target_id}'] = idmap_by_type[itype_key][target_id]
                except KeyError:
                    unmapped_ids.append(f'{itype_key}|{target_id}')

    with open(update_root / 'outputs/cache/target_result.json', 'w') as f:
        json.dump(result_dict, f, indent=2)
    
    with open(args.unmapped_log_file, 'w') as f:
        print(f'Dumping {len(unmapped_ids)} / {total_counts} unmapped IDs')
        f.write('\n'.join(unmapped_ids))
    # del idmap_by_type

    # print('Fetching protein description...')
    # if not prot_to_desc_json.exists() or OVERWRITE_SOURCE_MAPPING:
    #     target_fields = ['UniProtKB-ID', 'Gene_Name', 'Gene_ORFName', 'NCBI_TaxID']
    #     uprot2info_dict = fetch_protein_gene_desc(idmapping_source_fpath,
    #                                             prot_to_desc_json,
    #                                             target_fields=target_fields,
    #                                             chunksize=1000000)

    # if not prot_to_desc_tab.exists() or OVERWRITE_SOURCE_MAPPING:
    #     print('Re-write protein description into table...')
    #     with open(prot_to_desc_json, 'r') as f:
    #         uprot2info_dict = json.load(f)

    #     with open(prot_to_desc_tab, 'w') as f:
    #         f.write('\t'.join(['uprot']+target_fields) + '\n')
    #         for uprot, desc in uprot2info_dict.items():
    #             line = '\t'.join([desc[field] for field in target_fields])
    #             f.write(f'{uprot}\t{line}\n')

