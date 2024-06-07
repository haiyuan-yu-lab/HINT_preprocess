#!/usr/bin/env python
#DEPENDENCIES: 
#
#EMAIL_LOG:
#EMAIL_ERR:

'''Parse raw SIFTS files into PDB<->UniProt residue ranges'''

import sys, subprocess, glob, os, re, gzip
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.abspath('..'))

from config import *

from mjm_tools.mjm_parsers import zip_res_range

# sifts_output_dir = os.path.join(OUTPUT_DIR, 'sifts/parsed_files')
if not os.path.exists(SIFTS_PARSED_DIR):
    os.makedirs(SIFTS_PARSED_DIR, exist_ok=True)
output_file = os.path.join(SIFTS_MAPPING_FILE)
# sifts_data_path = '/home/resources/sifts/data/*/*.gz'
# output_file = '/home/resources/sifts/parsed_files/pdbresiduemapping.txt'
with open(PDB_TO_RUN) as f:
    pdb_list = f.read().splitlines()
sifts_files = ['{}/{}/{}.xml.gz'.format(SIFTS_DATA_DIR, pid.lower()[1:-1], pid.lower()) for pid in pdb_list]

#----------------------------------------------------------------------------
#parse what's already been calculated (don't repeat these to save time. To recalculate all, delete output file):
previously_calculated = defaultdict(list)
if os.path.exists(output_file):
    for l in open(output_file):
        previously_calculated[l.split('\t')[0]].append(l)   #lists of entries keyed on PDB

#----------------------------------------------------------------------------
# sifts_files = sorted(glob.glob(sifts_data_path), key=lambda f: os.path.basename(f))

output = open(output_file, 'w')
# YL: add taxonomy ID
output.write('\t'.join(['PDB', 'Chain', 'UniProt', 'TaxID', 'MappableResInPDBChainOnUniprotBasis', 'MappableResInPDBChainOnPDBBasis', 
                        'AllResInPDBChainOnPDBBasis', 'UniProtRes', 'UniProtResSecondaryStructure']) + '\n')


newly_calculated = 0
for s in sifts_files:
    pdb = os.path.basename(s)[:4].upper()
    
    #if it's already calculated, don't do it again
    if pdb in previously_calculated:
        for l in previously_calculated[pdb]:
            output.write(l)
        continue
    
    newly_calculated += 1
    # data = gzip.open(s).read()
    with gzip.open(s, 'rt') as f_xml:
        data = f_xml.read()
        residues = re.findall('\<residue dbSource.*?\>.*?\</residue\>', data, flags=re.DOTALL)
    
    chain2res = defaultdict(lambda: {'uniprot': '', 'taxa': '', 'upresnum': [], 'pdbresnum': [], 'allpdbres': [], 'upres': '', 'ss': ''})
    conflict_chains = set() #keep track of chains that have more than one or overlapping UniProt annotations	

    for r in residues:
    
        try:
            pdb, pdbresnum, chain = re.findall('\<crossRefDb dbSource="PDB" dbCoordSys="PDBresnum" dbAccessionId="(.*?)" dbResNum="(.*?)" dbResName=".*?" dbChainId="(.*?)"/\>', r)[0]
            chain2res[(pdb, chain)]['allpdbres'].append(pdbresnum)
        except IndexError:
            print(chain)
            continue
        
        
        #go no further if there is no match to uniprot residue for this pdb residue entry (examples: Conflict, Engineered mutation, Not_Observed) 
        if '<residueDetail dbSource="PDBe" property="Annotation">' in r:
            continue
        
        try:
            uniprot, upresnum, upres = re.findall('\<crossRefDb dbSource="UniProt" dbCoordSys="UniProt" dbAccessionId="(.*?)" dbResNum="(.*?)" dbResName="(.*?)"/>', r)[0]
        except IndexError:
            continue
        try:
            # dbSource="NCBI" dbCoordSys="UniProt" dbAccessionId="9606" dbResNum="20" dbResName="E"
            taxa = re.findall('\<crossRefDb dbSource="NCBI" dbCoordSys="UniProt" dbAccessionId="([0-9]+)" dbResNum=".*?" dbResName=".*?"/>', r)[0]
        except IndexError:
            taxa = "-"

        try:
            ss = re.findall('\<residueDetail dbSource="PDBe" property="codeSecondaryStructure">(.*?)\</residueDetail\>', r)[0]
        except IndexError:
            ss = "-"
        
        if chain2res[(pdb, chain)]['uniprot'] != '' and chain2res[(pdb, chain)]['uniprot'] != uniprot:
            conflict_chains.add((pdb,chain))	#keep track of chains that have more than one or overlapping UniProt annotations	

        chain2res[(pdb, chain)]['uniprot'] = uniprot
        chain2res[(pdb, chain)]['taxa'] = taxa
        chain2res[(pdb, chain)]['upresnum'].append(upresnum)
        chain2res[(pdb, chain)]['pdbresnum'].append(pdbresnum)
        chain2res[(pdb, chain)]['upres'] += upres
        chain2res[(pdb, chain)]['ss'] += ss
        
        
    for pdbchain in sorted(chain2res.keys()):

        if pdbchain in conflict_chains:
            continue
        
        pdb, chain = pdbchain

        upres = zip_res_range(chain2res[pdbchain]['upresnum'])
        pdbres = zip_res_range(chain2res[pdbchain]['pdbresnum'])
        allpdbres = zip_res_range(chain2res[pdbchain]['allpdbres'])
        
        if upres == '[]': continue
    
        output.write('\t'.join([pdb.upper(), chain, chain2res[pdbchain]['uniprot'], chain2res[pdbchain]['taxa'], 
                                upres, pdbres, allpdbres, chain2res[pdbchain]['upres'], chain2res[pdbchain]['ss']]) + '\n')
            

output.close()

print('Finished parsing SIFTS mappings for %s PDB structures (%s newly calculated).' %(len(sifts_files), newly_calculated))
