#!/usr/bin/env python
#DEPENDENCIES: pdb_03.py sifts_01.py
#
#EMAIL_LOG:
#EMAIL_ERR:

'''Roughly estimates the inter-chain distance for all PDB chain pairs to help decide which pairs it is worth running the full interface residue calculation on.'''

"""
Input: SIFTS_MAPPING_FILE ('{SIFTS_PARSED_DIR}/pdbresiduemapping.txt')
Output: {IRES_PARSED_DIR}/'bounding_box_pdb_chain_dists.txt'
"""

import os, sys, time
import gzip
import numpy as np
from collections import defaultdict
from itertools import combinations

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.abspath('..'))

from config import *

from mjm_tools.mjm_parsers import parse_dictionary_list
from mjm_tools.mjm_tools import time_elapsed
start_time = time.time()

#input files
# sifts_file = '/home/resources/sifts/parsed_files/pdbresiduemapping.txt'
#output files
if not os.path.exists(IRES_PARSED_DIR):
    os.makedirs(IRES_PARSED_DIR, exist_ok=True)

pdb_chain_dists_file = os.path.join(IRES_PARSED_DIR, 'bounding_box_pdb_chain_dists.txt')


#----------------------------------------------------------------------------
#parse what's already been calculated (don't repeat these to save time. To recalculate all, delete output file):

# start constructing chain combinations list so that ALL previously calculated results
# actually get re-written in the final output loop
# NOTE: Don't do this because the only entries that would be removed by not doing this
#       are obsolete PDB IDs. Defining the pdb2chaincombinations based on the up to date
#       sifts data ensures we don't retain obsolete entries. We WILL lose some records
#       between runs, but it will keep the resource up to date and clean.
#pdb2chaincombinations = defaultdict(set)


previously_calculated = {}
if os.path.exists(pdb_chain_dists_file):
    for l in open(pdb_chain_dists_file):
        uniA, uniB, pdb, cA, cB = tuple(l.strip().split('\t')[:5])
        #pdb2chaincombinations[pdb].add((uniA, uniB, pdb, cA, cB))
        previously_calculated[tuple(l.strip().split('\t')[:5])] = l   #keyed on (UniProtA, UniProtB, PDB, ChainA, ChainB)

#----------------------------------------------------------------------------
#get all chain combinations
sifts_data = parse_dictionary_list(SIFTS_MAPPING_FILE)

pdb2chainuniprot = defaultdict(set)
for e in sifts_data:
    pdb2chainuniprot[e['PDB']].add((e['UniProt'], e['Chain']))

pdb2chaincombinations = defaultdict(set)
total_combos = 0
for pdb in pdb2chainuniprot:
    chainpairs = combinations(pdb2chainuniprot[pdb], 2)
    for cp in [sorted(c) for c in chainpairs]:  #sort by uniprot, then chain
        (uniprotA, chainA), (uniprotB, chainB) = cp
        pdb2chaincombinations[pdb].add((uniprotA, uniprotB, pdb, chainA, chainB))
        total_combos += 1

#-----------------------------------------------------------------------------
#calculate bounding box based chain distances when necessary
pdbchains2dist = dict()
parsed_pdbs = 0
new_pairs1 = 0
print('Calculating chain distances in %s pdbs...' %(len(pdb2chaincombinations)))
for pdb in pdb2chainuniprot.keys():
    # Check to see if there are any combos in this chain that haven't already been calculated
    # otherwise, skip this PDB to save time
    new = len(pdb2chaincombinations[pdb].difference(previously_calculated))
    if(new == 0):
        continue
    
    # Try to open the file (may fail because of irregular pdb SIFTS entries)
    # pdb_fh = open_pdb(pdb, verbose=False)
    # if(pdb_fh is None):
    #     continue
    try:
        with gzip.open(os.path.join(PDB_DATA_DIR, '{}/pdb{}.ent.gz'.format(pdb.lower()[1:-1], pdb.lower())), 'rt') as pdb_fh:
            pdb_lines = pdb_fh.read().splitlines()
    except FileNotFoundError:
        continue
    
    # Counters for final output stats
    parsed_pdbs += 1
    new_pairs1 += new
    
    # Keep track of coordinates per chain
    chain2coors = defaultdict(lambda: {"xs":[], "ys":[], "zs":[]})
    for l in pdb_lines:
        # Skip any rows that aren't atom of hetatm (technically probably only care about atom but keep
        # hetatm to be safe)
        if(not(l[:6] in ["ATOM  ", "HETATM"])):
            continue
        chain = l[21:22]
        chain2coors[chain]["xs"].append(float(l[30:38]))
        chain2coors[chain]["ys"].append(float(l[38:46]))
        chain2coors[chain]["zs"].append(float(l[46:54]))
    
    # Generate min / max coordaintes for bounding box of each chain
    chains = chain2coors.keys()
    chain2bounding_box = dict()
    for chain in chains:
        mins = np.array([min(chain2coors[chain]["xs"]), min(chain2coors[chain]["ys"]), min(chain2coors[chain]["zs"])])
        maxes = np.array([max(chain2coors[chain]["xs"]), max(chain2coors[chain]["ys"]), max(chain2coors[chain]["zs"])])
        chain2bounding_box[chain] = (mins, maxes)
    
    # Iterate over all chain combos to calculate distance
    for chain1 in chains:
        mins1 = chain2bounding_box[chain1][0]
        maxes1 = chain2bounding_box[chain1][1]
        
        for chain2 in chains:
            mins2 = chain2bounding_box[chain2][0]
            maxes2 = chain2bounding_box[chain2][1]
            
            # Calculate linear distance between each coordinate
            dists = np.array([abs(mins1 - mins2), abs(mins1 - maxes2), abs(maxes1 - mins2), abs(maxes1 - maxes2)]).min(axis=0)
            
            # Zero out the distance for any coordinates where the bounding boxes overlap
            dists[(maxes2 >= mins1) & (mins2 <= maxes1)] = 0
            
            # Calculate rough distance between bounding boxes
            dist = np.sqrt(sum(dists**2))
            
            # Save this distance
            pdbchains2dist[(pdb, chain1, chain2)] = dist



#-----------------------------------------------------------------------------
# Write final output
output = open(pdb_chain_dists_file, "w")
output.write('\t'.join(['UniProtA','UniProtB','PDB','ChainA', 'ChainB', 'Distance']) + '\n')

new_pairs2 = 0
print('Writting final output for %s chain combinations...' %(total_combos))
calculated_pairs = set(pdbchains2dist.keys())
for pdb in pdb2chaincombinations.keys():
    for c in pdb2chaincombinations[pdb]:
        uniprotA, uniprotB, pdb, chainA, chainB = c
        
        if chainA==chainB:  #account for SIFTS behavior when 2 UniProts mapped to same chain 
            continue
        
        #if it's already calculated, don't do it again
        if c in previously_calculated: #(UniProtA, UniProtB, PDB, ChainA, ChainB)
            output.write(previously_calculated[c])
            continue
        
        new_pairs2 += 1
        if((pdb, chainA, chainB) in calculated_pairs):
            output.write('\t'.join([uniprotA, uniprotB, pdb, chainA, chainB, str(pdbchains2dist[(pdb, chainA, chainB)])])+'\n')
        else:
            # I'm not sure how to handle this. One option would be to assume the chains are close enough that we should test them
            # Another option would be to just skip the pair. I THINK the only way this should happen is if the PDB chain pair comes
            # from an invalid or irregular pdb file (in which case we shouldn't be able to calculate ires on them anyway)
            continue
            #output.write('\t'.join([uniprotA, uniprotB, pdb, chainA, chainB, -1])+'\n')

output.close()


print('Finished calculating interchain distances for %s chain combinations from SIFTS '
      'among %s PDBs. %s (or %s) newly calculated. Time elapsed: %s ' %(total_combos, parsed_pdbs, new_pairs1, new_pairs2, time_elapsed(start_time)))
