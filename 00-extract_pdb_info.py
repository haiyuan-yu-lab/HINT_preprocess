import glob, os, sys
import copy
import gzip
import tarfile
from tempfile import mkdtemp
from shutil import rmtree

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.abspath('..'))

# from config import *

"""
Summarize PDB information
Input:
    Path to PDB raw data & PDB-bundle raw data
        - `pdb_data_dir`, (default: /share/yu/resources/local_resources/pdb/data/data/pdb)
        - `pdb_bundle_data_dir`, (default: /share/yu/resources/local_resources/pdb/data/data/pdb_like)
    Output directory: `output_root`
Output:
    PDB meta information saved in $output_root/{pdb_info.txt, pdblike_info.txt}
"""

def read_pdb_format(pdb_fpath):
    if not os.path.exists(pdb_fpath):
        raise FileNotFoundError
    
    open_func = gzip.open if str(pdb_fpath).endswith('gz') else open

    categories = [('date', ''), ('pmid', ''), ('chains', []), ('exper', ''), ('resol', ''), ('nummod', ''), ('taxids', [])]
    pdb_data_dict = dict(categories)
    with open_func(pdb_fpath, 'rt') as f:
        for line in f:
            line = line.strip('\n')
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
    import argparse
    parser = argparse.ArgumentParser(description='Process PDB files')
    parser.add_argument('--output_root', help='Root directory for output files')
    parser.add_argument('--pdb_data_dir', default='/share/yu/resources/local_resources/pdb/data/data/pdb', 
                        help='Directory containing PDB files')
    parser.add_argument('--pdb_bundle_data_dir', default='/share/yu/resources/local_resources/pdb/data/data/pdb_like', 
                        help='Directory containing PDB-like files')
    args = parser.parse_args()

    output_root = args.output_root
    # parsed_dir = f'{output_root}/pdb/parsed_files'
    # if not os.path.exists(parsed_dir):
    #     os.makedirs(parsed_dir, exist_ok=True)
    
    RUN_PDB = True
    RUN_PDB_BUNDLE = True
    pdb_data_dir = args.pdb_data_dir

    if RUN_PDB:
        output_file = f'{output_root}/pdb_info.txt'
        print('Extracting PDB meta data...')

        categories = [('date', ''), ('pmid', ''), ('chains', []), ('exper', ''), ('resol', ''), ('nummod', ''), ('taxids', [])]
        
        all_pdb_fpaths = glob.glob(f'{pdb_data_dir}/*/*.gz')
        print(f'Found {len(all_pdb_fpaths)} PDB files to process.')
        # output_file = open('/home/resources/pdb/parsed_files/pdb_info.txt', 'w')
        with open(output_file, 'w') as fo:
            header = '\t'.join(['pdb',] + [c[0] for c in categories])
            fo.write(header + '\n')
            for pdb_fpath in all_pdb_fpaths:
                pdb_id = os.path.splitext(os.path.basename(pdb_fpath))[0][3:7].upper()
                try:
                    pdb_data = read_pdb_format(pdb_fpath)
                    # data[pdb_id] = pdb_data
                except FileNotFoundError:
                    print('FileNotFoundError: {}'.format(pdb_fpath))
                    continue
            
            # for pdb in sorted(data.keys()):
                tokens = [pdb_id,]
                for field, _ in categories:
                    token = ''
                    if type(pdb_data[field]) == type(set()) or type(pdb_data[field]) == type([]):
                        if len(pdb_data[field])==0: token = ''
                        else: token = ';'.join(pdb_data[field])
                    elif type(pdb_data[field]) == type(''):
                        token = pdb_data[field]

                    tokens.append(token)
                fo.write('\t'.join(tokens)+'\n')

        # fo.close()
        print(f'Finish parsing pdb info fields')
    
    if RUN_PDB_BUNDLE:
        output_file = f'{output_root}/pdblike_info.txt'
        print('Extracting PDB-like meta data...')
        # categories = [('date', ''), ('pmid', '')]
        header = ['pdb', 'date', 'pmid']
        # data = {}

        with open(output_file, 'w') as fo:
            fo.write('\t'.join(header) + '\n')
            all_pdb_fpaths = glob.glob(f'{args.pdb_bundle_data_dir}/*/*/*.tar.gz')
            print(f'Found {len(all_pdb_fpaths)} PDB-like files to process.')

            for n, pdblike_fpath in enumerate(all_pdb_fpaths):
                pdb_id = os.path.splitext(os.path.basename(pdblike_fpath))[0].split('-')[0].upper()
                try:
                    pdb_dict = read_pdb_bundle(pdblike_fpath)
                    pdb_dict['pdb'] = pdb_id
                    fo.write('\t'.join([pdb_dict[k].strip() for k in header]) + '\n')

                except FileNotFoundError:
                    print('FileNotFoundError: {}'.format(pdblike_fpath))
                    continue
        print(f'Finish parsing pdb-like info fields')

