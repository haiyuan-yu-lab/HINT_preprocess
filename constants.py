# Update options

DATA_CONFIGS = {
    'runBio': True,
    'runMINT': False, # MINT Downloads are no longer working
    'runIREF': True,
    'runDIP': True,
    'runIntAct': True,
    'runHPRD': True,
    'runMIPS': True,
    'runPDB': True
}

IRES_FILE = '/home/yl986/data/IRES/parsed_files/20240616/ires_perpdb_alltax.txt'
PDB_INFO_FILE = '/home/yl986/data/pdb/20240616/parsed_files/pdb_info.txt'
IRES_PDB_BUNDLE_FILE = '/home/yl986/data/IRES/parsed_files/20240616/ires_perpdb_alltax_pdblike.txt'
PDB_BUNDLE_INFO_FILE = '/home/yl986/data/pdb/20240616/parsed_files/pdblike_info.txt'

# URLs
# bioGridURL = "http://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ALL-LATEST.mitab.zip"
BIOGRID_URL = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.mitab.zip"
INTACT_URL = "ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip"

# manually merged taxa in initial parsing --- more taxa IDs to merge in later steps
REDUNDANT_TAXA = {
#     '4932': '559292',  # Cerv
    '559292': '4932',
    '562': '83333',    # E. coli
    '574521': '83333',   # E. coli
    '284812': '4896'
}

# Global ID dictionary for FTP data - maps different ID types to their UNIPROT DB KEY
# values are valid ID types from uniprot idmapping.dat file
ACCEPTED_IDS = {'refseq': "RefSeq", 
                'biogrid': "BioGRID", 
                'chembl target': "ChEMBL", 
                'complex': "ComplexPortal", 
                'ddbj/embl/genbank': "EMBL",
                'dip': "DIP",
                'ensembl': "Ensembl", 
                'ensemblgenomes': "EnsemblGenome",
                'entrezgene/locuslink': "GeneID",
                'entrez gene/locuslink': "GeneID",
                'reactome': "Reactome",
                'uniprot/swiss-prot': "UniProtKB-Swiss-Prot", 
                'uniprotkb': "UniProtKB",
                'wwpdb': "PDB"}

# Global ID dictionary for REST API - maps different ID types to their UNIPROT DB KEY
accepted_api_IDs = {'RefSeq': "RefSeq_Protein",
                    'refseq': "RefSeq_Protein", 
                    'biogrid': "BioGRID", 
                    'chembl target': "ChEMBL", 
                    'complex': "ComplexPortal", 
                    'ddbj/embl/genbank': "EMBL-GenBank-DDBJ",
                    'dip': "DIP",
                    'ensembl': "Ensembl", 
                    'ensemblgenomes': "Ensembl_Genomes",
                    'entrezgene/locuslink': "GeneID",
                    'entrez gene/locuslink': "GeneID",
                    'reactome': "Reactome",
                    'uniprot/swiss-prot': "UniProtKB-Swiss-Prot", 
                    'uniprotkb': "UniProtKB",
                    'wwpdb': "PDB"}