#!/usr/bin/env python
#DEPENDENCIES: download pdb files
#
#EMAIL_LOG:
#EMAIL_ERR:

import glob, os, sys
import copy
import gzip

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.abspath('..'))

from config import *

""""
Summarize PDB information
Input:
    PDB structure files for IDs specified in `PDB_TO_RUN`
Output:
    PDB meta information saved in `OUTPUT_DIR`/pdb/parsed_files/pdb_info.txt
"""

def read_pdb_format(pdb_fpath):
    if not os.path.exists(pdb_fpath):
        raise FileNotFoundError
    
    if str(pdb_fpath).endswith('gz'):
        with gzip.open(pdb_fpath, 'rt') as f:
            file_data = f.read().splitlines()
    else:
        with open(pdb_fpath) as f:
            file_data = f.read().splitlines()
    categories = [('date', ''), ('pmid', ''), ('chains', []), ('exper', ''), ('resol', ''), ('nummod', ''), ('taxids', [])]
    pdb_data_dict = dict(categories)

    for line in file_data:
        
        if 'SOURCE' == line[0:6] and 'ORGANISM_TAXID:' == line[11:26]:
            pdb_data_dict['taxids'].append(line[27:].split()[0].replace(';', ''))
        
        if 'HEADER' == line[0:6]:
            pdb_data_dict['date'] = line[50:59].split()[0]
        
        if 'JRNL        PMID' == line[0:16]:
            pdb_data_dict['pmid'] = line.strip().split()[-1]
        
        if 'REMARK   2 RESOLUTION.' == line[0:22]:
            pdb_data_dict['resol'] = line[22:].strip().replace('ANGSTROMS.', '').replace('NOT APPLICABLE.', 'n/a').strip()
        
        if 'EXPDTA' == line[0:6]:
            pdb_data_dict['exper'] = line[6:].strip().replace('X-RAY DIFFRACTION', 'XRAY').replace('SOLUTION NMR', 'NMR')
            
        if 'NUMMDL' == line[0:6]:
            pdb_data_dict['nummod'] = line[6:].strip()
            
        if 'COMPND' == line[0:6] and 'CHAIN:' == line[11:17]:
            pdb_data_dict['chains']+= line[17:].strip().replace(';', '').split(', ')
    
    return pdb_data_dict


if __name__ == '__main__':    
    # YL modified for testing
    output_root = OUTPUT_DIR
    parsed_dir = f'{output_root}/pdb/parsed_files'
    if not os.path.exists(parsed_dir):
        os.makedirs(parsed_dir, exist_ok=True)
    output_file = f'{parsed_dir}/pdb_info.txt'

    # all_pdb_files = sorted(glob.glob('/home/resources/pdb/data/*/*.gz'))

    with open(PDB_TO_RUN) as f:
        pdb_list = f.read().splitlines()
    all_pdb_files = ['{}/{}/pdb{}.ent.gz'.format(PDB_DATA_DIR, pid.lower()[1:-1], pid.lower()) for pid in pdb_list]

    categories = [('date', ''), ('pmid', ''), ('chains', []), ('exper', ''), ('resol', ''), ('nummod', ''), ('taxids', [])]
    data = {}

    for pdb_fpath in all_pdb_files:
        pdb_id = os.path.splitext(os.path.basename(pdb_fpath))[0][3:7].upper()
        try:
            pdb_data = read_pdb_format(pdb_fpath)
            data[pdb_id] = pdb_data
        except FileNotFoundError:
            print('FileNotFoundError: {}'.format(pdb_fpath))
            continue
        
    # output_file = open('/home/resources/pdb/parsed_files/pdb_info.txt', 'w')
    with open(output_file, 'w') as fo:
        header = '\t'.join(['pdb',] + [c[0] for c in categories])
        fo.write(header + '\n')

        for pdb in sorted(data.keys()):
            tokens = [pdb,]
            for field, _ in categories:
                token = ''
                if type(data[pdb][field]) == type(set()) or type(data[pdb][field]) == type([]):
                    if len(data[pdb][field])==0: token = ''
                    else: token = ';'.join(data[pdb][field])
                elif type(data[pdb][field]) == type(''):
                    token = data[pdb][field]
                
                tokens.append(token)
            fo.write('\t'.join(tokens)+'\n')

    # fo.close()
    print('parsed pdb info fields for %s pdb entries' %(len(all_pdb_files)))

