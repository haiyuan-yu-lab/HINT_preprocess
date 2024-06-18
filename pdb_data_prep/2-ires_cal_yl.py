'''
Batch calculate all interface residues in all SIFTS uniprot-tagged PDB structures, each pair of chains
Inputs: 
    - SIFTS_MAPPING_FILE
    - pdb_chain_dists_file (default: {IRES_PARSED_DIR}/'bounding_box_pdb_chain_dists.txt')
Outputs: 
    - IRES per PDB: (default: {IRES_PARSED_DIR}/'ires_perpdb_alltax.txt')
    - negative records in IRES calculation: {IRES_NEG_DIR}/NOires_perpdb_alltax.txt
'''

import sys, subprocess, os, time
import argparse
from collections import defaultdict
from itertools import combinations
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.abspath('..'))

from config import *

from mjm_tools.mjm_parsers import zip_res_range, unzip_res_range, parse_dictionary_list
from mjm_tools.mjm_tools import time_elapsed
from irescalc import ires_calc_fn


if __name__ == '__main__':
    ## I/O arguments at command line
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Calculates interface residues')
    parser.add_argument('-fdist', default='bounding_box_pdb_chain_dists.txt', help='file of bounding box distance between pdb chains')
    parser.add_argument('-output', default='ires_perpdb_alltax.txt', help='Output file (IRES per PDB among all taxonomy)', required=False)
    parser.add_argument('-fneg', default='NOires_perpdb_alltax.txt', help='Negative result file (IRES per PDB among all taxonomy)', required=False)

    args = parser.parse_args()

    start_time = time.time()

    # input_files
    sifts_file = SIFTS_MAPPING_FILE
    pdb_chain_dists_file = os.path.join(IRES_PARSED_DIR, args.fdist)

    # output_files
    if not os.path.exists(IRES_NEG_DIR):
        os.makedirs(IRES_NEG_DIR, exist_ok=True)

    negative_results_file = os.path.join(IRES_NEG_DIR, args.fneg)
    ires_output_fpath = os.path.join(IRES_PARSED_DIR, args.output)
    NEG_MODE = 'a'  # append
    OUTPUT_MODE = 'w'  # re-write

    # Distance threshold for which pairs get interface residues calculated
    close_thresh = DIST_THRES

    #----------------------------------------------------------------------------
    #parse what's already been calculated (don't repeat these to save time. To recalculate all, delete output file):
    previously_calculated = {}
    if os.path.exists(ires_output_fpath) and not os.path.getsize(ires_output_fpath) == 0:
        ires_prev = pd.read_csv(ires_output_fpath, sep='\t')
        if len(ires_prev):
            ires_prev['entry_key'] = ires_prev.apply(lambda x: (x['UniProtA'], x['UniProtB'], x['PDB'], x['ChainA'], x['ChainB']), axis=1)
            previously_calculated = set(ires_prev['entry_key'])
            ires_prev.drop('entry_key', axis=1).to_csv(ires_output_fpath, sep='\t', index=False, mode=OUTPUT_MODE)
            OUTPUT_MODE = 'a'
        # for l in open(ires_output_fpath):
        #     previously_calculated[tuple(l.strip().split('\t')[:5])] = l   #keyed on (UniProtA, UniProtB, PDB, ChainA, ChainB)
    

    #parse what couldn't be calculated before (don't try them again)
    previously_negative = set()
    if os.path.exists(negative_results_file) and not os.path.getsize(negative_results_file) == 0:
        # for l in open(negative_results_file):
        #     previously_negative.add(tuple(l.strip().split('\t')[:5]))  #(UniProtA, UniProtB, PDB, ChainA, ChainB)
        neg_prev = pd.read_csv(negative_results_file, sep='\t', names=['UniProtA','UniProtB','PDB','ChainA', 'ChainB'])
        if len(neg_prev):
            neg_prev['entry_key'] = neg_prev.apply(lambda x: (x['UniProtA'], x['UniProtB'], x['PDB'], x['ChainA'], x['ChainB']), axis=1)
            previously_negative = set(neg_prev['entry_key'])

    #----------------------------------------------------------------------------
    #get all chain combinations
    sifts_data = parse_dictionary_list(sifts_file)

    pdb2chainuniprot = defaultdict(set)
    pdbchainuniprot2sifts = dict()

    for entry in sifts_data:
        pdb2chainuniprot[entry['PDB']].add((entry['UniProt'], entry['Chain']))
        pdbchainuniprot2sifts[(entry['PDB'], entry['Chain'], entry['UniProt'])] = entry

    # chaincombinations = set()
    # for pdb in pdb2chainuniprot:
    #     chainpairs = combinations(pdb2chainuniprot[pdb], 2)
    #     for cp in [sorted(c) for c in chainpairs]:  #sort by uniprot, then chain
    #         (uniprotA, chainA), (uniprotB, chainB) = cp
    #         chaincombinations.add((uniprotA, uniprotB, pdb, chainA, chainB))
    #-----------------------------------------------------------------------------
    #parse chain distances
    df_dist = pd.read_csv(pdb_chain_dists_file, sep='\t')
    df_dist['entry_key'] = df_dist.apply(lambda x: (x['UniProtA'], x['UniProtB'], x['PDB'], x['ChainA'], x['ChainB']), axis=1)
    print('Calculating ires for %s chain combinations...' %(len(df_dist)))
    prev_result_indice = df_dist['entry_key'].isin(previously_calculated)
    prev_neg_indice = df_dist['entry_key'].isin(previously_negative)
    print('Skipping previously calculated results ({}) and previous negative entries ({})'.format(prev_result_indice.sum(), prev_neg_indice.sum()))
    df_dist = df_dist[~prev_result_indice & ~prev_neg_indice].reset_index(drop=True)
    print('Calculating {} new entries...'.format(df_dist.shape[0]))

    # pdbchains2dist = defaultdict(lambda: 100)   # Default distance set absurdly high if a distance has not been calculated
    # with open(pdb_chain_dists_file, "r") as f:
    #     for l in f:
    #         if("UniProt" in l): # Skip header
    #             continue
    #         pdbchains2dist[tuple(list(l.split("\t"))[:5])] = float(l.split("\t")[5])

    #-----------------------------------------------------------------------------
    #calculate ires when necessary

    # negative = open(negative_results_file, NEG_MODE)
    # output = open(ires_output_fpath, OUTPUT_MODE)
    with open(ires_output_fpath, OUTPUT_MODE) as fo:
        if OUTPUT_MODE == 'w':
            fo.write('\t'.join(['UniProtA','UniProtB','PDB','ChainA', 'ChainB', 'TaxIDA', 'TaxIDB', 'NumIresA', 'NumIresB', 'UniProtIresA', 'UniProtIresB', 
                                'PDBIresA', 'PDBIresB', 'UniProtAllResA', 'UniProtAllResB', 'PDBAllResA', 'PDBAllResB']) + '\n')

        newly_calculated = 0
        
        # for c in sorted(list(chaincombinations)):
        for i, record in df_dist.iterrows():
            if (i + 1) % 1000 == 0:
                print('{} / {} processed'.format(i+1, len(df_dist)))
            # uniprotA, uniprotB, pdb, chainA, chainB = c

            # Skip chains that are too far apart to possibly have an interface
            # if pdbchains2dist[(uniprotA, uniprotB, pdb, chainA, chainB)] > close_thresh:
            #     continue
            if record['Distance'] > close_thresh:
                continue
            
            if record['ChainA'] == record['ChainB']: # account for SIFTS behavior when 2 UniProts mapped to same chain 
                continue
            
            pdb = record['PDB']
            chainA = record['ChainA']
            chainB = record['ChainB']
            uniprotA = record['UniProtA']
            uniprotB = record['UniProtB']
            #if this is a pdb-like file or the PDB file for some other reason doesn't exist skip it
            if not os.path.exists(os.path.join(PDB_DATA_DIR, "{0}/pdb{1}.ent.gz".format(pdb[1:3].lower(), pdb.lower()))):
                continue
            try:
                #---------------------------------------------------
                #print uniprotA, uniprotB, pdb, chainA, chainB
                #chainA sifts residue mapping dictionary	
                pdbres = unzip_res_range(pdbchainuniprot2sifts[(pdb, chainA, uniprotA)]['MappableResInPDBChainOnPDBBasis'])
                uniprotres = unzip_res_range(pdbchainuniprot2sifts[(pdb, chainA, uniprotA)]['MappableResInPDBChainOnUniprotBasis'])
                
                chain_a_res2uniprot = dict(zip(pdbres, uniprotres))
                
                #chainB sifts residue mapping dictionary	
                # chainBres2uniprot = {}
                pdbres = unzip_res_range(pdbchainuniprot2sifts[(pdb, chainB, uniprotB)]['MappableResInPDBChainOnPDBBasis'])
                uniprotres = unzip_res_range(pdbchainuniprot2sifts[(pdb, chainB, uniprotB)]['MappableResInPDBChainOnUniprotBasis'])
                
                chain_b_res2uniprot = dict(zip(pdbres, uniprotres))
                
                #---------------------------------------------------
                
                # get interface residues
                # out, err = subprocess.Popen(['irescalc.py', pdb, '-c1', chainA, '-c2', chainB], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                try:
                    pdb_fpath = os.path.join(PDB_DATA_DIR, '{}/pdb{}.ent.gz'.format(pdb.lower()[1:-1], pdb.lower()))
                    PDBIresA, PDBIresB = ires_calc_fn(pdb_fpath, chainA, chainB, u_sasa=15, d_sasa=1)
                    # iresA, iresB, _ = out.split('\n')
                except ValueError as e:
                    print(e)
                    if "One or both" in str(e):
                        with open(negative_results_file, NEG_MODE) as f_neg:
                            f_neg.write('\t'.join(record['entry_key']) + '\n')
                            NEG_MODE = 'a'
                except Exception as e:
                    print(e)
                    continue
                
                # PDBIresA = iresA.split(',')
                # PDBIresB = iresB.split(',')
                
                UniProtIresA = [chain_a_res2uniprot[r] for r in PDBIresA if r in chain_a_res2uniprot]
                UniProtIresB = [chain_b_res2uniprot[r] for r in PDBIresB if r in chain_b_res2uniprot]
                
                if len(UniProtIresA) == 0 or len(UniProtIresB) == 0:
                    with open(negative_results_file, NEG_MODE) as f_neg:
                        f_neg.write('\t'.join(record['entry_key']) + '\n')
                        NEG_MODE = 'a'
                    continue
                
                nA = len(PDBIresA)
                nB = len(PDBIresB)
                PDBIresA = zip_res_range(PDBIresA)
                PDBIresB = zip_res_range(PDBIresB)
                UniProtIresA = zip_res_range(UniProtIresA)
                UniProtIresB = zip_res_range(UniProtIresB)
                
                UniProtAllResA = pdbchainuniprot2sifts[(pdb, chainA, uniprotA)]['MappableResInPDBChainOnUniprotBasis']
                UniProtAllResB = pdbchainuniprot2sifts[(pdb, chainB, uniprotB)]['MappableResInPDBChainOnUniprotBasis']
                
                PDBAllResA = pdbchainuniprot2sifts[(pdb, chainA, uniprotA)]['AllResInPDBChainOnPDBBasis']
                PDBAllResB = pdbchainuniprot2sifts[(pdb, chainB, uniprotB)]['AllResInPDBChainOnPDBBasis']
                
                taxidA = pdbchainuniprot2sifts[(pdb, chainA, uniprotA)]['TaxID']
                taxidB = pdbchainuniprot2sifts[(pdb, chainB, uniprotB)]['TaxID']

                fo.write('\t'.join([uniprotA, uniprotB, pdb, chainA, chainB, taxidA, taxidB, str(nA), str(nB), UniProtIresA, UniProtIresB, 
                                    PDBIresA, PDBIresB, UniProtAllResA, UniProtAllResB, PDBAllResA, PDBAllResB])+'\n')
                newly_calculated += 1
            
            except Exception as e:
                print(e)
                continue

    # negative.close()
    # output.close()

    print('Finished calculating interface residues for %s chain combinations from SIFTS '
        '(%s newly calculated with found ires). Time elapsed: %s ' %(len(df_dist), newly_calculated, time_elapsed(start_time)))

        
