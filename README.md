# HINT_preprocess
Data pre-processing pipeline for HINT

## Downloads

* Download up-to-date UniProt reference files from FTP site (https://ftp.uniprot.org/pub/databases/uniprot/current_release/)

  ```
  sh download.sh
  ```
Downloaded files will be saved at `$SOURCE_DIR` specified in the script

```
# SOURCE_DIR structure
.
|-- uniparc
|-- knowledgebase
    |-- idmapping
    |   |-- idmapping.dat.gz
    |-- complete
    |   |-- uniprot_sprot.fasta.gz  
    |   |-- uniprot_sprot_varsplic.fasta.gz  
    |   |-- uniprot_trembl.fasta.gz
    |-- docs
```

## Data parsing and processing

### 1. Parse UniProt reference files
* Run `parse_source_data.py` to process reference files downloaded from UniProt FTP site
    ```
    python parse_source_data.py $SOURCE_DIR/knowledgebase
    ```
    The following files will be processed:
    * FASTA files (`uniprot_sprot.fasta.gz`, `uniprot_sprot_varsplic.fasta.gz`, `uniprot_trembl.fasta.gz`) - to extract protein meta information
    * species file `docs/speclist.txt` to extract species taxonomy information
    * secondary-to-primary accession mapping file `docs/sec_ac.txt`

### 2. Parse and aggregate protein interaction dataset from different sources

Run `prepare_dataset.py` to prepare protein interaction data from sources of interest. The following steps will be processed by running the script:
1. Collect datasets from sources of interest. Active datasets will be downloaded from the source and inactive dataset will be copied from our self-maintained cache directory. Raw data files will be saved at `$UPDATE_DIR/data/parseTargets`

The following data sources are included (**Last update: 2024.4.15**)
* Active:
  * [BioGRID](https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.mitab.zip)
  * [IntAct](https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip)
  * PDB: copy to `static_datasets` from self-maintained PDB interaction data directory
* Inactive:
  * DIP (`dip20170205.txt`)
  * iRef (`All.mitab.03022013.txt`)
  * HPRD (`BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt`)
  * MIPS (`mppi` xml format)
  ```
  |-- static_datasets
      |-- dip20170205.txt
      |-- All.mitab.03022013.txt
      |-- mppi
      |-- BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt
      |-- ires_perpdb_alltax.txt
  ```
2. Parse raw data and generate initial `raw_interactions.txt` file saved at `$UPDATE_DIR/outputs/`

3. Revise parsed raw_interaction data and fill in UniProt IDs if available in source.Output files will be saved at `$UPDATE_DIR/outputs/cache/`. The following files are generated in this step:
  * `raw_interactions_filled_partial.txt` revised raw interaction files with UniProt IDs filled when available in source
  * `mapping_targets_by_type.json` IDs remain to be mapped to UniProt IDs orgainized by ID types. Supported ID types can be found in `constants.py`

### 3. Generate source-to-uniprot ID mapping

Run `create_idmapping.py` to parse `idmapping.dat.gz` from UniProt FTP site and generate ID mapping dictionary from source IDs to UniProt IDs

**Inputs**
* `idmapping.dat.gz`
* `mapping_targets_by_type.json`

**Outputs**
* `target_type_to_uprot.json` dictionary of ID mapping organized by ID type (if an ID is mapped to multiple UniProt IDs, all UniProt IDs will be kept & concatenated by `'|'`)
* `prot_gene_info.tsv` descriptions for each UniProt ID columns: `(uprot | UniProtKB-ID | Gene_Name | Gene_ORFName | NCBI_TaxID)`
* `target_result.json` mapping of target IDs to UniProt when available  
Format: `"ID_TYPE|SOURCE_ID": "UNIPROT_ID1(|UNIPROT_ID2|UNIPROT_ID3...)"` (Example: '"DIP|DIP-17064N": "Q9TW27"')

### 4. Proceed to remaining HINT data curation pipeline

Run codes in jupyter notebook: `4-HINT_data_curation.ipynb`
