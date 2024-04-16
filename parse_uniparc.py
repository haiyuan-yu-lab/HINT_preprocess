import string
import os, sys
from pathlib import Path
from glob import glob
import gzip
import numpy as np
import pandas as pd
import re

import json

from collections import defaultdict
from lxml import etree
import xml.etree.ElementTree as ET


def parse_entry(elem, ns):
    accession = elem.find('ns:accession', ns).text
    dbReferences = elem.findall('ns:dbReference', ns)
    try:
        seq_length = elem.find('ns:sequence', ns).get('length')
    except AttributeError:
        seq_length = ''
    acc_info_lst = []
    for dbRef in dbReferences:
        id_type = dbRef.get('type')
        if id_type in ["UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL"]:
            uprot = dbRef.get('id')
            is_reviewed = (id_type.find('TrEMBL') == -1)
            status = dbRef.get('active')
            try:
                taxa = dbRef.find('ns:property[@type="NCBI_taxonomy_id"]', ns).get('value')
            except AttributeError:
                taxa = ''
            try:
                protein_name = dbRef.find('.//ns:property[@type="protein_name"]', ns).get('value')
            except AttributeError:
                protein_name = ''
            try:
                gene_name = dbRef.find('.//ns:property[@type="gene_name"]', ns).get('value')
            except AttributeError:
                gene_name = ''
            cur_dict = {
                'accession': accession,
                'uniprot': uprot,
                'is_reviewed': is_reviewed,
                'status': status,
                'taxa': taxa,
                'protein_name': protein_name,
                'gene_name': gene_name,
                'seq_length': seq_length
            }
            acc_info_lst.append(cur_dict)
    return acc_info_lst


def parse_uniparc_xml(xml_path, output_fpath, inactive_fname='inactive.txt'):
    # Namespace
    ns = {'ns': 'http://uniprot.org/uniparc'}
    header = ['accession', 'uniprot', 'is_reviewed', 'status', 
              'taxa', 'protein_name', 'gene_name', 'seq_length', 'source_file']
    # if not overwrite and output_fpath.exists():
    #     print(str(output_fpath), 'already exists!')
    #     return
        
    # i = 0
    if not output_fpath.exists():
        with open(output_fpath, 'w') as f_out:
            f_out.write('\t'.join(header) + '\n')
    if not (uniparc_root / f'cache/{inactive_fname}').exists():
        with open(uniparc_root / f'cache/{inactive_fname}', 'w') as f_inac:
            f_inac.write('\t'.join(header) + '\n')

    f_inac = open(uniparc_root / f'cache/{inactive_fname}', 'a')
    if xml_path.name.endswith('gz'):
        f = gzip.open(xml_path, 'rb')
    else:
        f = open(xml_path)
    
    source_file = xml_path.name.split('.')[0]

    with open(output_fpath, 'a') as f_out:
        # f_out.write('\t'.join(header) + '\n')
        for event, elem in ET.iterparse(f, events=('start', 'end')):
            if event == 'end' and elem.tag.endswith('entry'):
                info_all = parse_entry(elem, ns)
                for tmp_dict in info_all:
                    tmp_dict['source_file'] = source_file
                    f_out.write('\t'.join([str(tmp_dict[k]) for k in header]) + '\n')
                    if tmp_dict['status'].upper() == 'N':
                        f_inac.write('\t'.join([str(tmp_dict[k]) for k in header]) + '\n')
                # Clear the element to free memory
                elem.clear()
                
    f.close()
    f_inac.close()


def extract_target_taxa(input_fpath, species_name2id, output_root, taxa_keyword='taxa', chunksize=1e6):
    """
    Extract target taxa from parsed uniparc data
    """
    if isinstance(input_fpath, str):
        input_fpath = Path(input_fpath)
    name = input_fpath.stem
    chunks = pd.read_csv(input_fpath, sep='\t', chunksize=chunksize, dtype=str)
    species_cache = set()
    output_header = True
    mode = 'w'
    for df in chunks:
        for species, taxa_lst in species_name2id.items():
            indices = df[taxa_keyword].isin(taxa_lst)
            species_tag = re.sub(' ', '', species.title()) 
            if indices.any():
                if species_tag in species_cache:
                    output_header = False
                    mode = 'a'
                df[indices].to_csv(output_root / f'{species_tag}_{name}.txt', index=False, sep='\t', mode=mode, header=output_header)
                species_cache.add(species_tag)


if __name__ == '__main__':
    uniparc_root = Path(sys.argv[1])
    uniparc_fname = sys.argv[2]
    output_root = Path(sys.argv[3])

    inactive_fname = 'inactive.txt'
    overwrite=True
    """
    Parse UniParc data
    """
    if not (uniparc_root / 'cache').exists():
        (uniparc_root / 'cache').mkdir(parents=True, exist_ok=True)

    uparc_output_fpath = uniparc_root / f'cache/uniparc_info_all.txt'

    xml_path = uniparc_root / uniparc_fname
    prefix = uniparc_fname.split('.')[0]
    uparc_output_fpath = uniparc_root / f'cache/{prefix}.txt'
    print(f'Processing {prefix}...')
    parse_uniparc_xml(xml_path, uparc_output_fpath, inactive_fname=inactive_fname)

    """
    Extract selected species from UniParc meta data 
    (uncomment to run the following when finished parsing uniparc raw file)
    """
    # with open(output_root / 'select_species2id.json') as f:
    #     select_species2id = json.load(f)
        
    # print('Extract selected species from UniParc meta data...')
    # output_root = uniparc_root / 'cache/taxa'
    # if not output_root.exists():
    #     output_root.mkdir(parents=True, exist_ok=True)

    # extract_target_taxa(uniparc_root / f'cache/{inactive_fname}', 
    #                     select_species2id, 
    #                     output_root)
    
    print('Done!')
