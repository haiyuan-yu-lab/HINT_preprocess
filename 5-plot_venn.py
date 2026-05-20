import os, sys
from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

"""
Generate Venn diagrams between new version and last active version (for human and yeast)
"""


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Compare HINT data versions')
    parser.add_argument('--hint_new_dir', type=str, help='Path to the update directory (parent of HINT format data)')
    parser.add_argument('--hint_prev_dir', type=str, help='Path to the previous directory (parent of HINT format data)')
    parser.add_argument('--timetag_new', type=str, help='Timestamp for the new version')
    parser.add_argument('--timetag_old', type=str, help='Timestamp for the old version')

    args = parser.parse_args()
    hint_prev_root = Path(args.hint_prev_dir)
    hint_new_root = Path(args.hint_new_dir)
    # ------- Path configurations --------
    time_tag_old = args.timetag_old
    time_tag_new = args.timetag_new
    # archive_root = data_root / f'old_versions/{time_stamp_old}/HINT_format/taxa/'

    # hint_output_root = output_root / 'HINT_format'
    # new_data_root = hint_output_root / 'taxa'

    fig_dir = hint_new_root / 'figures'
    if not fig_dir.exists():
        fig_dir.mkdir()
    
    # suffix of target file to process
    target_suffix = ['binary_all',
                    'binary_hq',
                    'both_all',
                    'both_hq',
                    'cocomp_all',
                    'cocomp_hq',
                    'htb_hq',
                    'htc_hq',
                    'lcb_hq',
                    'lcc_hq']

    target_species = ['Homo sapiens', 'Saccharomyces cerevisiae']
    target_species = [s.lower() for s in target_species]

    for species in target_species:
        species_tag = re.sub(' ', '', species.title()) 
        
        for suffix in target_suffix:
            try:
                ppi_new = pd.read_csv(hint_new_root / f'HINT_format/taxa/{species_tag}/{species_tag}_{suffix}.txt', sep='\t')
                ppi_new['ppi'] = ppi_new.apply(lambda x: ':'.join(sorted([x['Uniprot_A'], x['Uniprot_B']])), axis=1)
            except FileNotFoundError as e:
                print(e)
                continue
            try:
                ppi_old = pd.read_csv(hint_prev_root / f'HINT_format/taxa/{species_tag}/{species_tag}_{suffix}.txt', sep='\t')
                ppi_old['ppi'] = ppi_old.apply(lambda x: ':'.join(sorted([x['Uniprot_A'], x['Uniprot_B']])), axis=1)
            except FileNotFoundError as e:
                print(e)
                continue

            # generate venn diagrams
            counts_old = ppi_old['ppi'].nunique()
            counts_new = ppi_new['ppi'].nunique()
            plt.figure()
            plt.title(species)
            g = venn2([set(ppi_old['ppi']), set(ppi_new['ppi'])], (f'{time_tag_old}-{suffix}\n({counts_old:,})', 
                                                                   f'{time_tag_new}-{suffix}\n({counts_new:,})'))
            plt.savefig(fig_dir / f'{species_tag}_{suffix}_{time_tag_old}_{time_tag_new}.pdf', bbox_inches='tight')