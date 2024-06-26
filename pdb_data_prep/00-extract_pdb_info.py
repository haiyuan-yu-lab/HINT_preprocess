#!/usr/bin/env python
#DEPENDENCIES: download pdb files
#
#EMAIL_LOG:
#EMAIL_ERR:

import glob, os, sys
import copy
import gzip
import tarfile
from tempfile import mkdtemp
from shutil import rmtree

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
        
        if line.startswith('ATOM') or line.startswith('HETATM'):
            break
    
    return pdb_data_dict


def read_pdb_bundle(pdb_fpath):
    # Create a scratch space for extracting intermediary files
    scratchDir = mkdtemp()

    # Extract PDB-lik data to scratch directory
    # fh = tarfile.open(pdb_like_f)
    categories = [('date', ''), ('pmid', '')]
    pdb_data_dict = dict(categories)
    pmid_retrieved = False
    
    with tarfile.open(pdb_fpath) as fh:
        fh.extractall(scratchDir)
        pdb_bundles = [x for x in fh.getnames() if not "chain" in x]

    for bundle in pdb_bundles:
        if pmid_retrieved:
            break
        for line in open("{0}/{1}".format(scratchDir, bundle), "r"):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                break
            if 'HEADER' == line[0:6]:
                pdb_data_dict['date'] = line[50:60].split()[0]
            
            if 'JRNL        PMID' == line[0:16]:
                pdb_data_dict['pmid'] = line.strip().split()[-1]
                pmid_retrieved = True
                break
                
    rmtree(scratchDir)
    
    return pdb_data_dict

if __name__ == '__main__':    
    output_root = OUTPUT_DIR
    parsed_dir = f'{output_root}/pdb/parsed_files'
    if not os.path.exists(parsed_dir):
        os.makedirs(parsed_dir, exist_ok=True)
    
    RUN_PDB = False
    RUN_PDB_BUNDLE = True

    if RUN_PDB:
        output_file = f'{parsed_dir}/pdb_info.txt'

        # all_pdb_files = sorted(glob.glob('/home/resources/pdb/data/*/*.gz'))
        if not os.path.exists(PDB_TO_RUN):
            all_pdb_files = sorted(glob.glob(f'{PDB_DATA_DIR}/*/*.gz'))
            pdb_list = [os.path.basename(s).split('.')[0][3:] for s in all_pdb_files]
            with open(PDB_TO_RUN, 'w') as f:
                f.write('\n'.join(pdb_list))

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
    
    if RUN_PDB_BUNDLE:
        output_file = f'{parsed_dir}/pdblike_info.txt'
        with open(PDB_BUNDLE_LIST) as f:
            pdb_list = f.read().splitlines()
        all_pdb_files = ['{}/{}/{}-pdb-bundle.tar.gz'.format(PDB_BUNDLE_DATA_DIR, pid.lower()[1:-1], pid.lower()) for pid in sorted(pdb_list)]

        # categories = [('date', ''), ('pmid', '')]
        header = ['pdb', 'date', 'pmid']
        # data = {}

        with open(output_file, 'w') as fo:
            fo.write('\t'.join(header) + '\n')
            for pdblike_fpath in all_pdb_files:
                pdb_id = os.path.splitext(os.path.basename(pdblike_fpath))[0].split('-')[0].upper()
                try:
                    pdb_dict = read_pdb_bundle(pdblike_fpath)
                    pdb_dict['pdb'] = pdb_id
                    fo.write('\t'.join([pdb_dict[k].strip() for k in header]) + '\n')

                except FileNotFoundError:
                    print('FileNotFoundError: {}'.format(pdblike_fpath))
                    continue
        print('parsed info fields for %s pdb-like entries' %(len(all_pdb_files)))

