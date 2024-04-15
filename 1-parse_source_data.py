import os, sys
import re
from pathlib import Path

from glob import glob
import gzip
import numpy as np
import pandas as pd

from collections import defaultdict

from Bio import SeqIO

"""
Usage: python parse_source_data.py /PATH/TO/UNIPROT/KNOWLEDGEBASE
"""

def species_parser(data_root, fname='speclist.txt', start_idx=59):
    """
    parse `speclist.txt` from UniProt FTP site
    """
    # skip the headers (set row index manually)
    taxa_list = []
    name_types = {'N': 'scientific_name', 'C': 'common_name', 'S': 'synonym'}
    with open(data_root / fname, 'r') as f:
        lines = f.read().splitlines()

    cur_dict = dict()
    for i, line in enumerate(lines[start_idx:]):
        if not line:
            break
        if line.startswith('======'):  # reach the end of species section
            break
        if line[0] != " ":
            if cur_dict:
                taxa_list.append(cur_dict)
            
            code, letter, taxa, desc = line[:6].strip(), line[6].strip(), line[7:16].strip(': '), line[16:].strip()
            cur_dict = {'code': code, 'kingdom': letter, 'taxa': taxa}
            desc_type, desc_info = desc.split('=')
        else:
            desc_type, desc_info = line[16:].strip().split('=')
        cur_dict[name_types[desc_type]] = desc_info

    taxa_list.append(cur_dict)
    df_species = pd.DataFrame.from_dict(taxa_list)

    return df_species


def sec_ac_parser(data_root, fname='sec_ac.txt', start_idx=31):
    """
    Parse `sec_acc.txt` (secondary to primary accession mapping file) from UniProt FTP site
    """
    if isinstance(data_root, str):
        data_root = Path(data_root)

    with open(data_root / fname, 'r') as f:
        lines = f.read().splitlines()

    if not (data_root / 'cache').exists():
        (data_root / 'cache').mkdir(parents=True)
        
    with open(data_root / 'cache/{}_parsed.txt'.format(fname.split('.')[0]), 'w') as f:
        f.write('\t'.join(['secondary', 'primary']) + '\n')
        for line in lines[start_idx:]:
            f.write(re.sub(r'\s+', '\t', line) + '\n')


def extract_info_fields(record):
    tag, pid, entry_name = record.id.split('|')
    prot_dict = {'tag': tag, 'UniProt': pid, 'name': entry_name}
    prot_dict['length'] = len(record.seq)
    desc = record.description.split(' ', 1)[1]  # remove record.id part
    os_search = re.search(r'OS=', desc)
    if os_search:
        prot_desc_end, os_start = os_search.span()
        prot_dict['description'] = desc[:prot_desc_end].strip()
    else:
        prot_dict['description'] = ''
    
    prot_dict['species'] = ''
    ox_search = re.search(r'OX=', desc)
    if ox_search:
        os_end, ox_start = ox_search.span()
        if os_search:
            prot_dict['species'] = desc[os_start:os_end].strip()
    try:
        prot_dict['taxa'] = re.search(r'OX=(.*?)\s', desc).group(1)
    except AttributeError:
        prot_dict['taxa'] = ''
    try:
        prot_dict['gene_name'] = re.search(r'GN=(.*?\S+)\s?', desc).group(1)
    except AttributeError:
        prot_dict['gene_name'] = ''
    
    return prot_dict


if __name__ == '__main__':
    source_dir = sys.argv[1]

    source_root = Path(source_dir)
    fasta_root = source_root / 'complete'
    doc_root = source_root / 'docs'
    prot_meta_root = fasta_root / 'meta'

    if not prot_meta_root.exists():
        prot_meta_root.mkdir(parents=True, exist_ok=True)
    
    """
    Extract protein meta information from fasta files:
    Inputs:
        - uniprot_sprot.fasta.gz
        - uniprot_sprot_varsplic.fasta.gz
        - uniprot_trembl.fasta.gz
    Outputs:
        - meta/sprot_meta.txt
        - meta/sprot_varsplic_meta.txt
        - meta/trembl_meta.txt
    """
    fields = ['tag', 'UniProt', 'name', 'length', 'species', 'taxa', 'gene_name', 'description']
    
    fasta_files = list(fasta_root.glob('*fasta.gz'))

    for fpath in fasta_files:
        fasta_name = fpath.name
        out_prefix = fasta_name.split('.')[0].split('_', 1)[-1]
        print(f'Processing {fasta_name}...')
        fa_records = SeqIO.parse(gzip.open(fpath, 'rt'), 'fasta')

        with open(prot_meta_root / f'{out_prefix}_meta.txt', 'w') as f:
            f.write('\t'.join(fields)+'\n')
            for record in fa_records:
                try:
                    prot_dict = extract_info_fields(record)
                    f.write('\t'.join(str(prot_dict[k]) for k in fields)+'\n')
                except:
                    print(record.description)
    
    doc_cache_root = doc_root / 'cache'
    if not doc_cache_root.exists():
        doc_cache_root.mkdir(parents=True, exist_ok=True)

    """
    Parse species file: 
    Input:
        $SOURCE_ROOT/docs/speclist.txt
    Output:
        $SOURCE_ROOT/docs/cache/species_parsed.txt
    """
    print('Parsing species file...')
    df_species = species_parser(doc_root, 'speclist.txt', start_idx=59)
    df_species.to_csv(doc_cache_root / 'species_parsed.txt', index=False, sep='\t')

    """
    Parse secondary-to-primary accession mapping file
    Input:
        $SOURCE_ROOT/docs/sec_ac.txt
    Output:
        $SOURCE_ROOT/docs/cache/sec_ac_parsed.txt
    """
    print('Parsing secondary-to-primary mapping file...')
    sec_ac_parser(doc_root, 'sec_ac.txt', start_idx=31)
