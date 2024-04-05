import os, sys
from pathlib import Path
import subprocess
import string
import urllib
# import urllib2
import datetime
import time
import shutil
from glob import glob
import gzip
import pandas as pd
from bs4 import BeautifulSoup 
from bs4 import BeautifulStoneSoup


sys.path.append(os.path.dirname(os.path.dirname(__file__)))
# sys.path.append(os.path.dirname(__file__))
# sys.path.append(os.path.abspath('..'))

from constants import *

# Wrapper function to download files given an URL, used to load datasets
def downloadFile(url, target_dir):
    os.system('wget -P {} -q {}'.format(target_dir, url))
    
def get_url_contents(url):
    request = urllib.request.urlopen(url)
    content = request.read()
    request.close()
    return content


def download_datasets(data_root):
    """
    PIPELINE STEP 1: Dataset downloading pipeline
    READS : data/static_datasets/*
    WRITES: data/parseTargets/*
        Notes made by Junke on 10/09/2023
        * DIP seems to have stopped updating since 20170205
        * Add 2 more static datasets since some databases are no longer updating
    """
    if isinstance(data_root, str):
        data_root = Path(data_root)

    tmp_path = data_root / 'tmp'
    parse_target_path = data_root / 'parseTargets'
    static_data_root = data_root / 'static_datasets'
    if not parse_target_path.exists():
        parse_target_path.mkdir(parents=True)
    if not static_data_root.exists():
        static_data_root.mkdir()

    # Remove the previous parseTargets folder
    if 'parseTargets' in os.listdir(data_root):
        shutil.rmtree(parse_target_path)
    if 'tmp' in os.listdir(data_root):
        shutil.rmtree(tmp_path)

    # Create the folders we will be working on and enter tmp
    tmp_path.mkdir(parents=True)
    parse_target_path.mkdir(parents=True)
    # os.chdir('./data/tmp')

    # BioGrid Download and unzip
    # bioGridVersion = get_url_contents(
    #     'http://webservice.thebiogrid.org/version?accesskey=d1e4afe96b70c4c4668b07d8b7635cae').read()
    # if bioGridVersion in ['3.4.131', '3.4.132', '3.4.133']:
    #     bioGridVersion = '3.4.134'
    downloadFile(BIOGRID_URL, str(parse_target_path))
    biogrid_zip = str(parse_target_path / BIOGRID_URL.split('/')[-1])
    subprocess.call(['unzip', '-d', str(parse_target_path), biogrid_zip])
    # biogrid_txt_file = glob('BIOGRID*.txt')[0]
    # shutil.move(biogrid_txt_file, str(parse_target_path / biogrid_txt_file))

    # Intact Download
    
    downloadFile(INTACT_URL, str(parse_target_path))
    intact_zip = str(parse_target_path / INTACT_URL.split('/')[-1])
    subprocess.call(['unzip', '-d', str(parse_target_path), intact_zip])
    # shutil.move("intact.txt", parse_target_path / 'intact.txt')

    # MINT Download
    # mintFolder = "ftp://mint.bio.uniroma2.it/pub/release/mitab26/current/"
    # mraw_interactions = [x.strip().split(' ')[-1] for x in get_url_contents(mintFolder).split('\n') if 'full-binary' in x][0]
    # downloadFile(mintFolder + mraw_interactions)
    # shutil.move(mraw_interactions, '../parseTargets/%s' % (mraw_interactions))

    # downloadFile("http://mint.bio.uniroma2.it/mitab/MINT_MiTab.txt")
    # shutil.move('MINT_MiTab.txt', '../parseTargets/%s' % ('MINT_MiTab.txt'))

    # # MIPS Download and unzip
    # mipsURL = "http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz"
    # downloadFile(mipsURL)
    # subprocess.call(['gzip', '-d', 'mppi.gz'])
    # shutil.move("mppi", '../parseTargets/mppi')

    # iRefWeb Download and unzip
    # iRefWebURL = "ftp://ftp.no.embnet.org/irefindex/data/archive/"
    # files = [x.split(' ')[-1]
    #          for x in get_url_contents(iRefWebURL).split('\n') if 'release' in x]
    # iRefWebRelease = sorted(
    #     files, key=lambda x: float(x.split('_')[-1]))[-1].strip()
    # iRefSubDir = iRefWebURL + iRefWebRelease + "/psi_mitab/MITAB2.6/"
    # iRefFile = [x.split(' ')[-1].strip()
    #             for x in get_url_contents(iRefSubDir).split('\n') if 'All' in x][0]
    # downloadFile(iRefSubDir + iRefFile)
    # subprocess.call(['unzip', iRefFile])
    # iRef = iRefFile.replace('.zip', '')
    # shutil.move(iRef, '../parseTargets/%s' % (iRef))

    # Move the static datasets (Which are placed there by a user, not automatically)
    # os.chdir('../..')
    staticFiles = os.listdir(static_data_root)
    for fname in staticFiles:
        shutil.copy(static_data_root / fname, parse_target_path / fname)

    # Now that we are done, remove the tmp folder
    shutil.rmtree(tmp_path)

    # Print reminder to admin about updating static_datasets
    print("=" * 20)
    print("REMEMBER TO MANUALLY REPLACE FILES IN data/static_datasets FOR:")
    print("""
        ================
        DIP database   =
        ================

        URL: http://dip.doe-mbi.ucla.edu/dip/Download.cgi?SM=3
        u: jfbeltran p: yulab

        FILE: dip<date>.txt
        + Tab separated (MITAB)

        ================
        Human Protein Reference Database
        ================

        URL: http://www.hprd.org/download
        SUBMIT AGREEMENT FORM
        UNZIP: tar zxvf HPRD_Release<n>_<date>.tar.gz
        FILE: HPRD_Release<n>_<date>/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt
        + Tab separated
        
        ================
        PDB
        ================
        
        copy parsed files from 
        /fs/cbsuhyfs1/storage/resources/ires/parsed_files/ires_perpdb_alltax.txt &
        /fs/cbsuhyfs1/storage/resources/ires/parsed_files/pdb_info.txt
        to static_datasets
        """)


def standardizeMITAB(fpath, sourceName, filehandler):
    """
    Reads a MITAB-formatted File and writes to a smaller 
    TAB-separated version of the data to filehandler
    """
    print("\nParsing MITAB file for {0} at {1}... ".format(sourceName, str(fpath)))

    f = open(fpath, 'r')
    header = f.readline().strip('#' + string.whitespace).split('\t')

    for line in f.readlines():
        # Increase linecount
#         print(ROW_COUNTER)
        ROW_COUNTER[0] += 1
        if ROW_COUNTER[0] % 50000 == 0:
            print("\n[%s] Working on ROW %s" % (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), ROW_COUNTER[0]))

        # Get rid of noise around lines
        row = line.strip('#' + string.whitespace).split('\t')

        # Grab the fields
        first_ID_field, second_ID_field, method_field, publication_field, first_taxa_field, second_taxa_field = ([(subrow.split(':')[0], ':'.join(
            subrow.split(':')[1:])) if subrow.count(':') > 1 else subrow.split(':') for subrow in row[i].split('|')] for i in (0, 1, 6, 8, 9, 10))

        # We do not accept publication serial numbers in DIP or imex format
        publications = [p[1] for p in publication_field if p[
            0] == 'pubmed' and 'DIP' not in p[1] and 'unassigned' not in p[1]]
        methods = [m[-1].split('(')[0].strip('"').split(':')[1] if ':' in m[-1].split('(')[0].strip('"') else m[-1].split('(')[0].strip('"') for m in method_field]
        methods = list(filter(lambda x: x != '-', methods))

        # When there is only one method listed, we assume it applies to all
        # publications
        if len(methods) == 1 and len(publications) > 1:
            methods = methods * len(publications)

        # When there is only one publication listed, we assume it uses all the
        # methods listed
        if len(publications) == 1 and len(methods) > 1:
            publications = publications * len(methods)

        # We do not accept publications without a matching methods
        if len(publications) != len(methods):
            if publications and methods:
                print('=' * 20)
                print("WARNING - Uneven indexes for:\nPUBLICATIONS:{0}\nMETHODS:{1}".format(publications, methods))
                print(publication_field, method_field)
                print('=' * 20)
        publication_methods = zip(publications, methods)

        # We do not accept inter-taxa interactions
        # taxid:83334(Escherichia coli O157:H7)
        taxas = set([taxa[-1].split('(')[0].strip('"')
                     for taxa in first_taxa_field + second_taxa_field])
        taxas = set(
            [taxa if taxa not in REDUNDANT_TAXA else REDUNDANT_TAXA[taxa] for taxa in taxas])
        if len(taxas) > 1:
            continue
        taxa = taxas.pop()

        for publicationID, methodID in publication_methods:
            for idA in first_ID_field:
                for idB in second_ID_field:
                    values = [
                        '|'.join(idA), '|'.join(idB), methodID, publicationID, taxa]
                    filehandler.write(
                        '\t'.join([sourceName, str(ROW_COUNTER[0])] + values) + '\n')
    f.close()
    

def standardize_datasets(
        update_root,
        runBio = True,
        runMINT = False, # MINT Downloads are no longer working
        runIREF = True,
        runDIP = True,
        runIntAct = True,
        runHPRD = True,
        runMIPS = True,
        runPDB = True):
    """
    PIPELINE STEP 2: Parse all datasets into one big TAB-separated file.
    READS : data/parseTargets/*
    WRITES: outputs/raw_interactions.txt
    """
    if isinstance(update_root, str):
        update_root = Path(update_root)
    data_root = update_root / 'data'
    parse_target_path = update_root / 'data/parseTargets'
    pdb_data_root = data_root / 'static_datasets'
    output_root = update_root / 'outputs'
    if not output_root.exists():
        output_root.mkdir(parents=True)
    # Initialize metadata for rows
    head = ["source", "row_number", "idA", "idB", "method", "publ", "taxa"]
    final = open(output_root / 'raw_interactions.txt', 'w')
    final.write('\t'.join(head) + '\n')

    # Identify the corresponding file or files for each dataset
    files = os.listdir(parse_target_path)
    print(files)
    # biogrid = "data/parseTargets/" + \
    #     [f for f in files if "BIOGRID" in f and 'txt' in f][0]
    # mint = "data/parseTargets/" + [f for f in files if "MINT" in f][0]
    # iref = "data/parseTargets/" + [f for f in files if "All.mitab" in f][0]
    # dip = "data/parseTargets/" + [f for f in files if f[:3] == 'dip'][0]
    # intact = "data/parseTargets/" + [f for f in files if f == 'intact.txt'][0]
    # hprd = "data/parseTargets/" + \
    #     [f for f in files if "BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt" == f][0]
    # mips = "data/parseTargets/" + [f for f in files if "mppi" == f][0]

    # CUSTOMIZABLE: Set flags for which datasets you would like to process
    # Parsing of each MITAB dataset
    if runDIP:
        try:
            dip = list(parse_target_path.glob('dip*txt'))[0]
            standardizeMITAB(dip, "DIP", final)
        except IndexError:
            print('DIP data not found')
        
    if runBio:
        try:
            biogrid = list(parse_target_path.glob('BIOGRID*txt'))[0]
            standardizeMITAB(biogrid, "BioGrid", final)
        except IndexError:
            print('BioGRID data not found')
            runBio = False
        
    if runMINT:
        try:
            mint = list(parse_target_path.glob('*MINT*'))[0]
            standardizeMITAB(mint, "MINT", final)
        except IndexError:
            print('MINT data not found')
        
    if runIREF:
        try:
            iref = list(parse_target_path.glob('All.mitab*txt'))[0]
            standardizeMITAB(iref, "iRef", final)
        except IndexError:
            print('iRef data not found')
            runIREF = False
        
    if runIntAct:
        try:
            intact = list(parse_target_path.glob('intact.txt'))[0]
            standardizeMITAB(intact, "IntAct", final)
        except IndexError:
            print('IntAct data not found')
            runIntAct = False
        
    # HPRD Specific:
    if runHPRD:
        experiment_codes = {'y2h': '0018', 'yeast 2-hybrid': '0018',
                                'in vivo': '0493', 'in vitro': '0492', 'vv': '0493', 'vt': '0492'}
        hprd = parse_target_path / 'BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt'
        if not hprd.exists():
            print('HPRD data not found')
        else:
            print("\nparsing HPRD ")
            # 9606 for human taxa
            f = open(hprd, 'r')
            for line in f.readlines():
                ROW_COUNTER[0] += 1
                if ROW_COUNTER[0] % 50000 == 0:
                    print("\n[%s] Working on ROW %s" % (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), ROW_COUNTER[0]))

                row = line.strip().split('\t')
                for experiment_type in row[6].split(';'):
                    for pubmedID in row[7].split(','):
                        idA, idB = '|'.join(["RefSeq", row[2]]), '|'.join(
                            ["RefSeq", row[5]])
                        readout = ['HPRD', str(ROW_COUNTER[0]), idA, idB, experiment_codes[
                            experiment_type.strip().lower()], pubmedID, "9606"]
                        final.write('\t'.join(readout) + "\n")

    # MIPS Specific:
    if runMIPS:
        mips = parse_target_path / 'mppi'
        if not mips.exists():
            print('MIPS data not found')
        else:
            print("\nparsing MIPS XML ")
            try:
                mppi = BeautifulSoup(open(mips).read(), features="html.parser")
                allInteractions = mppi.findAll('interaction')
            except:
                mppi = BeautifulStoneSoup(open(mips).read())
                allInteractions = mppi.findAll('interaction')

            for i, interaction in enumerate(allInteractions):
                ROW_COUNTER[0] += 1
                if ROW_COUNTER[0] % 50000 == 0:
                    print("\n[%s] Working on ROW %s" % (datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'), ROW_COUNTER[0]))

                try:
                    primaries = [r['id']
                                for r in interaction.findAll('primaryref')]
                    organisms = [o['ncbitaxid']
                                for o in interaction.findAll('organism')]
                    organisms = [
                        o if o not in REDUNDANT_TAXA else REDUNDANT_TAXA[o] for o in organisms]
                    idA, idB = '|'.join(["uniprotkb", primaries[2]]), '|'.join(["uniprotkb", primaries[3]])
                except KeyError:
                    continue

                if len(set(organisms)) > 1:
                    continue

                readout = ['MIPS', str(ROW_COUNTER[0]), idA, idB, primaries[1], primaries[0], organisms[0]]
                final.write('\t'.join(readout) + "\n")

    if runPDB:
        print("\nparsing PDB Logs ")
        if not ((pdb_data_root / 'ires_perpdb_alltax.txt').exists() and (pdb_data_root / 'pdb_info.txt').exists()):
            print('PDB data not available')
        else:
            interactions = open(pdb_data_root / 'ires_perpdb_alltax.txt', 'r')
            interactions.readline()
    #         publications = dict([(line.split('\t')[0], (line.split('\t')[2], line.split('\t')[4])) for line in open(
    #             '/home/resources/pdb/parsed_files/pdb_info.txt').read().strip().split('\n')[1:]])
            publications = dict([(line.split('\t')[0], (line.split('\t')[2], line.split('\t')[4])) \
                                for line in open(pdb_data_root / 'pdb_info.txt').read().strip().split('\n')[1:]])
            dumped = set()
            for line in interactions.readlines():
                ROW_COUNTER[0] += 1
                interaction = line.split('\t')
                idA, idB, pdb = interaction[:3]
                taxaA, taxaB = interaction[5:7]
                # Merge redundant taxa
                if(taxaA in REDUNDANT_TAXA):
                    taxaA = REDUNDANT_TAXA[taxaA]
                if(taxaB in REDUNDANT_TAXA):
                    taxaB = REDUNDANT_TAXA[taxaB]
                
                if pdb not in publications:
                    dumped.add(line)
                    continue
                pubID, experiment_type = publications[pdb]
                if pubID.strip() == '':
                    pubID = 'PDB_' + str(pdb)
                if taxaA != taxaB:
                    dumped.add(line)
                    continue
                readout = ['PDB', str(ROW_COUNTER[0]), 'uniprotkb|' + idA,
                        'uniprotkb|' + idB, '0114', pubID, taxaA]
                final.write('\t'.join(readout) + "\n")
            # MODIFIED BY SHAYNE FOR DEBUGGING PURPOSES
            #with open('/home/resources/interactomes/outputs/dumped_pdbs.txt', 'w') as f:
            with open(output_root / 'dumped_pdbs.txt', 'w') as f:
                f.write('\n'.join([str(x) for x in sorted(list(dumped))]))
    final.close()
    print("\nFinished Standardization of Datasets\n")


def replace_explicit_interactome(update_root):
    """
    PIPELINE STEP 3: Replace all entries from certain publication-method-taxa combos

    READS: 
      * `data/explicit_interactomes/*.txt`
      * `outputs/raw_interactions.txt`

    WRITES: `outputs/raw_interactions.txt`
    """
    if isinstance(update_root, str):
        update_root = Path(update_root)
    
    data_root = update_root / 'data'
    output_root = update_root / 'outputs'
    if not output_root.exists():
        output_root.mkdir(parents=True, exist_ok=True)

    print("Replacing poor interactions using explicit interactome")
    # explicit_files = glob('data/explicit_interactomes/*.txt')
    explicit_files = list((data_root / 'explicit_interactomes').glob('*txt'))

    # Load the whole interactions file into memory like savages
    raw_interactions = [line.split('\t') for line in open(output_root / 'raw_interactions.txt').read().strip().split('\n')]
    ROW_COUNTER[0] = int(raw_interactions[-1][1])

    for fpath in explicit_files:
        # We grab the pub, taxa, and method details from the filename
        filename = fpath.name
        publication, taxa, method = filename.split('.')[0].split('_')
        
        # Merge redundant taxa
        if(taxa in REDUNDANT_TAXA):
            taxa = REDUNDANT_TAXA[taxa]
            
        # Create the new lines
        newLines = [[filename, str(ROW_COUNTER[0] + i + 1), '\t'.join(['uniprotkb|' + uniprot for uniprot in ids.split('\t')]), method, publication, taxa]
                    for i, ids in enumerate(open(fpath).read().strip().split('\n'))]

        # Increase the counter
        ROW_COUNTER[0] += len(newLines)

        # Delete all redundant lines from the raw_interactions
        raw_interactions = [r for r in raw_interactions if r[3] != method and r[4] != publication and r[5] != taxa]

        # Add the new lines
        raw_interactions.extend(newLines)

    # Write the whole thing back to the file
    with open(output_root / 'raw_interactions.txt', 'w') as f:
        for line in raw_interactions:
            f.write('\t'.join(line) + '\n')  
    

if __name__ == '__main__':
    # Set up starting / output directories
    # home = os.getcwd()
    update_dir = "/home/yl986/data/HINT/update_2024"
    archive_data_dir = '/home/yl986/data/HINT/data/'
    wipe = False
    ROW_COUNTER = [0]

    """
    Step 0 - set up paths
    """
    # Recreate the directories from scratch
    if(wipe):
        os.system("rm -rf {0}".format(update_dir))

    # Create output directory if it does not already exist
    update_root = Path(update_dir)
    if not update_root.exists():
        update_root.mkdir(parents=True)
        # os.system("mkdir {0}".format(update_dir))
    data_root = update_root / 'data'
    if not data_root.exists():
        data_root.mkdir(parents=True)
    
    # Copy over necessary files that are not automatically created
    # Static Datasets
    os.system("cp {}/static_datasets/ {}/data/ -r".format(archive_data_dir, update_dir))

    # Bouncer
    os.system("cp {}/bouncer {}/data/ -r".format(archive_data_dir, update_dir))

    # Explicit Interactomes
    os.system("cp {}/explicit_interactomes/ {}/data/ -r".format(archive_data_dir, update_dir))

    # Taxon ID to Names
    os.system("cp {}/taxid2name.txt {}/data/".format(archive_data_dir, update_dir))

    # Evidence Codes
    os.system("cp {}/evidence_code*.txt {}/data/".format(archive_data_dir, update_dir))

    # Create any necessary file structure not handled by the pipeline itself
 
    # Outputs folder
    # os.system("mkdir {0}/outputs".format(output_dir))

    # # HINT_format folder
    # os.system("mkdir {0}/outputs/HINT_format".format(output_dir))

    # # taxa folder
    # os.system("mkdir {0}/outputs/HINT_format/taxa".format(output_dir))

    # # Caches folder
    # os.system("mkdir {0}/caches".format(output_dir))

    # os.chdir(output_dir)
    # #os.chdir('/home/resources/interactomes')
    
    """
    Step 1 - Download datasets
    """
    download_datasets(data_root)
    
    """
    Step 2 - Parse raw data and generate initial `raw_interactions.txt` file
    
    Reads in all of the datasets in data/parseTargets, standardizes the format, and saved the output to outputs/raw_interactions.txt
    This step also reads in information from our IRES and PDB resources to find PDB interactions.
    Need to confirm that both of these are updating properly.
    """
    standardize_datasets(update_dir, **DATA_CONFIGS)
    print('Number of rows parsed:', ROW_COUNTER[0])

    """
    Step 3 - Replace all entries from certain publication-method-taxa combos
    Uses some pre-generated (don't know how, where, when, or why) interactome sets for some organisms (data/explicit_interactomes)
    Adds all of the entries included in these files to `raw_interactions.txt` and drops any of the redundant rows.
    """
    replace_explicit_interactome(update_dir)
    # Dataset ready for curation
