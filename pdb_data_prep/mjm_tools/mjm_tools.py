"""Contains useful tools for interacting with UniProt, PDB, and other biological data sources. """

import os
from gzip import open as gzopen
# import urllib, urllib2
# from urllib2 import urlopen, HTTPError
import urllib.parse
import urllib.request
from urllib.request import urlopen, HTTPError
import unicodedata
from collections import defaultdict
from re import split as regexsplit
from re import sub, findall
import time
from subprocess import Popen, PIPE
import warnings
import copy

from itertools import product

from scipy.stats import norm
import numpy as np
import math

from config import *

from mjm_tools.mjm_parsers import unzip_res_range, zip_res_range
from mjm_tools.mjm_parsers import parse_dictionary_list

# command used to run nacccess (either absolute path to the binary, or simply the program name 'naccess' if program is on a system path)
# NOTE: several files used by naccess must be in the same directory as the naccess executable: vdw.radii, standard.data, accall, etc.
# AND naccess needs to have been installed (using install.scr) in the directory in which it currently resides (dependency paths are hard-coded into the naccess binary).
# NACCESS_PATH = 'naccess'

# command used to run freesasa (either absolute path to the binary, or simply the program name 'naccess' if program is on a system path)
# NOTE: 2020_03_03 we are considering this as an alternative to naccess because we have found some cases where naccess seems to
#       erroneously calculates a dSASA between structures even when there should be no dSASA (e.g. irescalc.py -c1 A -c2 F 3NSW)
# FREESASA_PATH = '/home/sdw95/bin/freesasa'

# command used to run ProFit (either absolute path to the binary, or simply the program name 'profit' if program is on a system path)
# PROFIT_PATH = 'profit'

# location of mjm_tools (also containing mjm_parsers.py, irescalc.py, srescalc.py, crescalc.py, etc.)
# TOOLS_PATH = '/home/mjm659/mjm_path/'

# absolute path to local mirrored location of PDB data
# PDB_DATA_DIR = '/home/resources/pdb/data/'

#absolute path to SIFTS residue mapping file
# SIFTS_MAPPING_FILE = '/home/resources/sifts/parsed_files/pdbresiduemapping.txt'

###############################################################################
################################# Statistics ##################################
###############################################################################


def std_err(p, n):
    '''Standard error of a binomial.
    p: fraction observed (or probability)
    n: number of trials (or denominator of fraction)'''
    
    p = float(p); n = float(n)
    return ((p*(1-p))/n)**0.5


def bino2test(x1, n1, x2, n2):
    '''To test whether the percentages in two samples are the same;
    i.e. to test wether the p's in two binomial distributions are the same
    x1: number of successes in sample one;
    n1: number of total trials in sample one;
    x2: number of successes in sample two;
    n2: number of total trials in sample two;
    Converted to Python from code by Haiyuan Yu, 10/31/2004'''

    x1 = float(x1); n1 = float(n1); x2 = float(x2); n2 = float(n2)

    p1 = float(x1) / float(n1)
    p2 = float(x2) / float(n2)
    ph = (x1 + x2) / (n1 + n2) #% ph: p-hat, the estimated sample p under the null hypothesis that p1 == p2;
    z = abs(p1 - p2) / (ph*(1-ph)*(1/n1 + 1/n2))**0.5
    p = 2*norm.cdf(-1*z)
    return p


def cl_effect_size(sample1, sample2):
    '''Compute the Common Language effect size to accompany Mann Whitney U-test
    (https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test#Common_language_effect_size) '''
    med1, med2 = np.median(sample1), np.median(sample2)
    all_pw_combinations = list(product(sample1, sample2))
    
    num_greater = 0.0
    
    for p in all_pw_combinations:
        if med1 > med2:
            if p[0] > p[1]: num_greater += 1
        elif med2 > med1:
            if p[1] > p[0]: num_greater += 1
        else:
            return None
    
    return num_greater/len(all_pw_combinations)

###############################################################################
################################## Plotting ###################################
###############################################################################

#colors from www.tinygorilla.com/Easter_eggs/PallateHex.html (now defunct?)

pastels = ['#6ECFF6', '#FFF79A', '#BC8DBF', '#F7977A', '#82CA9D', '#C4DF9B', '#A187BE', '#F6989D', '#8493CA', '#FDC68A', '#7BCDC8', '#A2D39C', '#F49AC2', '#7EA7D8', '#8882BE']
lights =  ['#A67C52', '#855FA8', '#F26C4F', '#7CC576', '#438CCA', '#5574B9', '#00BFF3', '#A763A8', '#605CA8', '#FFF467', '#ACD372', '#F06EA9', '#FBAF5C', '#1ABBB4', '#3BB878', '#F26D7D']
darkers = ['#32004B', '#005B7F', '#005826', '#7B0046', '#827B00', '#37302D', '#790000', '#7A0026', '#4B0049', '#002157', '#603913', '#0D004C', '#005952', '#7D4900', '#005E20', '#003663', '#406618']
darks =   ['#598527', '#630460', '#007236', '#440E62', '#9E0039', '#9E005D', '#004B80', '#003471', '#0076A3', '#A36209', '#754C24', '#00746B', '#1A7B30', '#1B1464', '#ABA000', '#534741', '#9E0B0F']
grays =   ['#C2C2C2', '#CCCCCC', '#262626', '#A0A0A0', '#707070', '#464646', '#B7B7B7', '#111111', '#ACACAC', '#898989', '#555555', '#E1E1E1', '#7D7D7D', '#626262', '#D7D7D7', '#363636', '#959595']

#http://www.mulinblog.com/a-color-palette-optimized-for-data-visualization/
mjm_colors = ['#5DA5DA', '#F15854', '#60BD68', '#B276B2', '#FAA43A', '#4D4D4D', '#F17CB0', '#DECF3F', '#B2912F']


all_colors = {'darker_violet': '#32004B', 'cmyk_green': '#00A651', '30%_gray': '#C2C2C2', 'yellow': '#FFF200', 'medium_warm_brown': '#8C6239', '25%_gray': '#CCCCCC', 'dark_pea_green': '#598527', 'light_warm_brown': '#A67C52', '90%_gray': '#262626', 'light_violet': '#855FA8', 'rgb_blue': '#0000FF', 'light_red': '#F26C4F', '45%_gray': '#A0A0A0', '65%_gray': '#707070', 'yellow_green': '#39B54A', 'cyan': '#00AEEF', 'dark_violet_magenta': '#630460', 'darker_cyan': '#005B7F', 'pastel_cyan': '#6ECFF6', '80%_gray': '#464646', 'rgb_green': '#00FF00', 'pea_green': '#8DC73F', 'darker_green': '#005826', 'dark_green': '#007236', 'light_yellow_green': '#7CC576', 'green_cyan': '#00A99D', 'pastel_yellow': '#FFF79A', 'dark_violet': '#440E62', 'rgb_red': '#FF0000', 'violet': '#662D91', 'dark_magenta_red': '#9E0039', 'pastel_violet_magenta': '#BC8DBF', 'darker_magenta': '#7B0046', 'blue': '#0054A6', 'blue_violet': '#2E3192', 'dark_magenta': '#9E005D', 'pastel_red': '#F7977A', '35%_gray': '#B7B7B7', 'darker_yellow': '#827B00', 'darker_cool_brown': '#37302D', 'red': '#ED1C24', 'light_cyan_blue': '#438CCA', 'medium_cool_brown': '#736357', 'pastel_green': '#82CA9D', 'pastel_pea_green': '#C4DF9B', 'light_blue': '#5574B9', 'pastel_violet': '#A187BE', 'light_cyan': '#00BFF3', 'darker_red': '#790000', '95%_gray': '#111111', 'pastel_magenta_red': '#F6989D', 'dark_cyan_blue': '#004B80', '40%_gray': '#ACACAC', 'light_violet_magenta': '#A763A8', 'yellow_orange': '#F7941D', 'light_blue_violet': '#605CA8', 'pastel_blue': '#8493CA', 'cmyk_magenta': '#EC008C', 'light_yellow': '#FFF467', 'dark_blue': '#003471', 'dark_cyan': '#0076A3', 'pastel_yellow_orange': '#FDC68A', 'darker_magenta_red': '#7A0026', 'light_pea_green': '#ACD372', '55%_gray': '#898989', 'darker_violet_magenta': '#4B0049', 'pastel_green_cyan': '#7BCDC8', 'light_magenta': '#F06EA9', 'darker_blue': '#002157', '75%_gray': '#555555', 'cmyk_red': '#ED1C24', 'pale_warm_brown': '#C69C6E', '15%_gray': '#E1E1E1', '60%_gray': '#7D7D7D', 'dark_yellow_orange': '#A36209', 'pastel_yellow_green': '#A2D39C', 'darker_warm_brown': '#603913', 'cyan_blue': '#0072BC', 'darker_blue_violet': '#0D004C', 'magenta_red': '#ED145B', 'light_yellow_orange': '#FBAF5C', 'darker_green_cyan': '#005952', 'pastel_magenta': '#F49AC2', '70%_gray': '#626262', 'dark_warm_brown': '#754C24', 'darker_yellow_orange': '#7D4900', 'violet_magenta': '#92278F', 'magenta': '#EC008C', '20%_gray': '#D7D7D7', 'cmyk_yellow': '#FFF200', 'dark_green_cyan': '#00746B', 'dark_yellow_green': '#1A7B30', 'pastel_cyan_blue': '#7EA7D8', 'darker_yellow_green': '#005E20', 'darker_cyan_blue': '#003663', 'dark_blue_violet': '#1B1464', 'rgb_cyan': '#00FFFF', 'dark_yellow': '#ABA000', 'pastel_blue_violet': '#8882BE', 'light_green_cyan': '#1ABBB4', '85%_gray': '#363636', 'dark_cool_brown': '#534741', 'light_green': '#3BB878', 'darker_pea_green': '#406618', 'dark_red': '#9E0B0F', 'light_magenta_red': '#F26D7D', 'rgb_magenta': '#FF00FF', 'green': '#00A651', '50%_gray': '#959595', 'cmyk_blue': '#2E3192', 'cmyk_cyan': '#00AEEF'}

def printc(s, c=31):
    '''Print string s with UNIX formatting (colors and styles here: http://misc.flogisoft.com/bash/tip_colors_and_formatting) '''
    print('\x1b[1;%sm' %(c) + s + '\x1b[0m')  #default is red


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def adjust_hex_brightness(hex_color, brightness_offset=50):
    """ Lightens or darkens a hex color """
    rgb_hex = [hex_color.lstrip('#')[x:x+2] for x in [0, 2, 4]]
    new_rgb_int = [int(hex_value, 16) + brightness_offset for hex_value in rgb_hex]
    new_rgb_int = [min([255, max([0, i])]) for i in new_rgb_int] # make sure new values are between 0 and 255
    return "#" + "".join([hex(i)[2:].zfill(2) for i in new_rgb_int])


###############################################################################
################################# GENERAL BIO #################################
###############################################################################

aa_dict = defaultdict(lambda: 'X', {   #default to ambiguous amino acid, X
    'ALANINE':'A', 'ARGININE':'R', 'ASPARAGINE':'N', 'ASPARTICACID':'D',
    'CYSTEINE':'C', 'GLUTAMICACID':'E', 'GLUTAMINE':'Q', 'GLYCINE':'G',
    'HISTIDINE':'H', 'ISOLEUCINE':'I', 'LEUCINE':'L', 'LYSINE':'K',
    'METHIONINE':'M', 'PHENYLALANINE':'F', 'PROLINE':'P', 'SERINE':'S',
    'THREONINE':'T', 'TRYPTOPHAN':'W', 'TYROSINE':'Y', 'VALINE':'V',
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
    'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I',
    'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
    'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
    'A':'A', 'R':'R', 'N':'N', 'D':'D', 'C':'C',
    'E':'E', 'Q':'Q', 'G':'G', 'H':'H', 'I':'I',
    'L':'L', 'K':'K', 'M':'M', 'F':'F', 'P':'P',
    'S':'S', 'T':'T', 'W':'W', 'Y':'Y', 'V':'V',
    })
    
codon_dict = defaultdict(lambda: 'X', {   #default to ambiguous amino acid, X
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"
    })


grantham = {
('A', 'A'): 0, ('A', 'C'): 195, ('A', 'D'): 126, ('A', 'E'): 107, ('A', 'F'): 113, ('A', 'G'): 60, ('A', 'H'): 86, ('A', 'I'): 94, ('A', 'K'): 106, ('A', 'L'): 96, ('A', 'M'): 84, ('A', 'N'): 111, ('A', 'P'): 27, ('A', 'Q'): 91, ('A', 'R'): 112, ('A', 'S'): 99, ('A', 'T'): 58, ('A', 'V'): 64, ('A', 'W'): 148, ('A', 'Y'): 112,
('C', 'A'): 195, ('C', 'C'): 0, ('C', 'D'): 154, ('C', 'E'): 170, ('C', 'F'): 205, ('C', 'G'): 159, ('C', 'H'): 174, ('C', 'I'): 198, ('C', 'K'): 202, ('C', 'L'): 198, ('C', 'M'): 196, ('C', 'N'): 139, ('C', 'P'): 169, ('C', 'Q'): 154, ('C', 'R'): 180, ('C', 'S'): 112, ('C', 'T'): 149, ('C', 'V'): 192, ('C', 'W'): 215, ('C', 'Y'): 194,
('D', 'A'): 126, ('D', 'C'): 154, ('D', 'D'): 0, ('D', 'E'): 45, ('D', 'F'): 177, ('D', 'G'): 94, ('D', 'H'): 81, ('D', 'I'): 168, ('D', 'K'): 101, ('D', 'L'): 172, ('D', 'M'): 160, ('D', 'N'): 23, ('D', 'P'): 108, ('D', 'Q'): 61, ('D', 'R'): 96, ('D', 'S'): 65, ('D', 'T'): 85, ('D', 'V'): 152, ('D', 'W'): 181, ('D', 'Y'): 160,
('E', 'A'): 107, ('E', 'C'): 170, ('E', 'D'): 45, ('E', 'E'): 0, ('E', 'F'): 140, ('E', 'G'): 98, ('E', 'H'): 40, ('E', 'I'): 134, ('E', 'K'): 56, ('E', 'L'): 138, ('E', 'M'): 126, ('E', 'N'): 42, ('E', 'P'): 93, ('E', 'Q'): 29, ('E', 'R'): 54, ('E', 'S'): 80, ('E', 'T'): 65, ('E', 'V'): 121, ('E', 'W'): 152, ('E', 'Y'): 122,
('F', 'A'): 113, ('F', 'C'): 205, ('F', 'D'): 177, ('F', 'E'): 140, ('F', 'F'): 0, ('F', 'G'): 153, ('F', 'H'): 100, ('F', 'I'): 21, ('F', 'K'): 102, ('F', 'L'): 22, ('F', 'M'): 28, ('F', 'N'): 158, ('F', 'P'): 114, ('F', 'Q'): 116, ('F', 'R'): 97, ('F', 'S'): 155, ('F', 'T'): 103, ('F', 'V'): 50, ('F', 'W'): 40, ('F', 'Y'): 22,
('G', 'A'): 60, ('G', 'C'): 159, ('G', 'D'): 94, ('G', 'E'): 98, ('G', 'F'): 153, ('G', 'G'): 0, ('G', 'H'): 98, ('G', 'I'): 135, ('G', 'K'): 127, ('G', 'L'): 138, ('G', 'M'): 127, ('G', 'N'): 80, ('G', 'P'): 42, ('G', 'Q'): 87, ('G', 'R'): 125, ('G', 'S'): 56, ('G', 'T'): 59, ('G', 'V'): 109, ('G', 'W'): 184, ('G', 'Y'): 147,
('H', 'A'): 86, ('H', 'C'): 174, ('H', 'D'): 81, ('H', 'E'): 40, ('H', 'F'): 100, ('H', 'G'): 98, ('H', 'H'): 0, ('H', 'I'): 94, ('H', 'K'): 32, ('H', 'L'): 99, ('H', 'M'): 87, ('H', 'N'): 68, ('H', 'P'): 77, ('H', 'Q'): 24, ('H', 'R'): 29, ('H', 'S'): 89, ('H', 'T'): 47, ('H', 'V'): 84, ('H', 'W'): 115, ('H', 'Y'): 83,
('I', 'A'): 94, ('I', 'C'): 198, ('I', 'D'): 168, ('I', 'E'): 134, ('I', 'F'): 21, ('I', 'G'): 135, ('I', 'H'): 94, ('I', 'I'): 0, ('I', 'K'): 102, ('I', 'L'): 5, ('I', 'M'): 10, ('I', 'N'): 149, ('I', 'P'): 95, ('I', 'Q'): 109, ('I', 'R'): 97, ('I', 'S'): 142, ('I', 'T'): 89, ('I', 'V'): 29, ('I', 'W'): 61, ('I', 'Y'): 33,
('K', 'A'): 106, ('K', 'C'): 202, ('K', 'D'): 101, ('K', 'E'): 56, ('K', 'F'): 102, ('K', 'G'): 127, ('K', 'H'): 32, ('K', 'I'): 102, ('K', 'K'): 0, ('K', 'L'): 107, ('K', 'M'): 95, ('K', 'N'): 94, ('K', 'P'): 103, ('K', 'Q'): 53, ('K', 'R'): 26, ('K', 'S'): 121, ('K', 'T'): 78, ('K', 'V'): 97, ('K', 'W'): 110, ('K', 'Y'): 85,
('L', 'A'): 96, ('L', 'C'): 198, ('L', 'D'): 172, ('L', 'E'): 138, ('L', 'F'): 22, ('L', 'G'): 138, ('L', 'H'): 99, ('L', 'I'): 5, ('L', 'K'): 107, ('L', 'L'): 0, ('L', 'M'): 15, ('L', 'N'): 153, ('L', 'P'): 98, ('L', 'Q'): 113, ('L', 'R'): 102, ('L', 'S'): 145, ('L', 'T'): 92, ('L', 'V'): 32, ('L', 'W'): 61, ('L', 'Y'): 36,
('M', 'A'): 84, ('M', 'C'): 196, ('M', 'D'): 160, ('M', 'E'): 126, ('M', 'F'): 28, ('M', 'G'): 127, ('M', 'H'): 87, ('M', 'I'): 10, ('M', 'K'): 95, ('M', 'L'): 15, ('M', 'M'): 0, ('M', 'N'): 142, ('M', 'P'): 87, ('M', 'Q'): 101, ('M', 'R'): 91, ('M', 'S'): 135, ('M', 'T'): 81, ('M', 'V'): 21, ('M', 'W'): 67, ('M', 'Y'): 36,
('N', 'A'): 111, ('N', 'C'): 139, ('N', 'D'): 23, ('N', 'E'): 42, ('N', 'F'): 158, ('N', 'G'): 80, ('N', 'H'): 68, ('N', 'I'): 149, ('N', 'K'): 94, ('N', 'L'): 153, ('N', 'M'): 142, ('N', 'N'): 0, ('N', 'P'): 91, ('N', 'Q'): 46, ('N', 'R'): 86, ('N', 'S'): 46, ('N', 'T'): 65, ('N', 'V'): 133, ('N', 'W'): 174, ('N', 'Y'): 143,
('P', 'A'): 27, ('P', 'C'): 169, ('P', 'D'): 108, ('P', 'E'): 93, ('P', 'F'): 114, ('P', 'G'): 42, ('P', 'H'): 77, ('P', 'I'): 95, ('P', 'K'): 103, ('P', 'L'): 98, ('P', 'M'): 87, ('P', 'N'): 91, ('P', 'P'): 0, ('P', 'Q'): 76, ('P', 'R'): 103, ('P', 'S'): 74, ('P', 'T'): 38, ('P', 'V'): 68, ('P', 'W'): 147, ('P', 'Y'): 110,
('Q', 'A'): 91, ('Q', 'C'): 154, ('Q', 'D'): 61, ('Q', 'E'): 29, ('Q', 'F'): 116, ('Q', 'G'): 87, ('Q', 'H'): 24, ('Q', 'I'): 109, ('Q', 'K'): 53, ('Q', 'L'): 113, ('Q', 'M'): 101, ('Q', 'N'): 46, ('Q', 'P'): 76, ('Q', 'Q'): 0, ('Q', 'R'): 43, ('Q', 'S'): 68, ('Q', 'T'): 42, ('Q', 'V'): 96, ('Q', 'W'): 130, ('Q', 'Y'): 99,
('R', 'A'): 112, ('R', 'C'): 180, ('R', 'D'): 96, ('R', 'E'): 54, ('R', 'F'): 97, ('R', 'G'): 125, ('R', 'H'): 29, ('R', 'I'): 97, ('R', 'K'): 26, ('R', 'L'): 102, ('R', 'M'): 91, ('R', 'N'): 86, ('R', 'P'): 103, ('R', 'Q'): 43, ('R', 'R'): 0, ('R', 'S'): 110, ('R', 'T'): 71, ('R', 'V'): 96, ('R', 'W'): 101, ('R', 'Y'): 77,
('S', 'A'): 99, ('S', 'C'): 112, ('S', 'D'): 65, ('S', 'E'): 80, ('S', 'F'): 155, ('S', 'G'): 56, ('S', 'H'): 89, ('S', 'I'): 142, ('S', 'K'): 121, ('S', 'L'): 145, ('S', 'M'): 135, ('S', 'N'): 46, ('S', 'P'): 74, ('S', 'Q'): 68, ('S', 'R'): 110, ('S', 'S'): 0, ('S', 'T'): 58, ('S', 'V'): 124, ('S', 'W'): 177, ('S', 'Y'): 144,
('T', 'A'): 58, ('T', 'C'): 149, ('T', 'D'): 85, ('T', 'E'): 65, ('T', 'F'): 103, ('T', 'G'): 59, ('T', 'H'): 47, ('T', 'I'): 89, ('T', 'K'): 78, ('T', 'L'): 92, ('T', 'M'): 81, ('T', 'N'): 65, ('T', 'P'): 38, ('T', 'Q'): 42, ('T', 'R'): 71, ('T', 'S'): 58, ('T', 'T'): 0, ('T', 'V'): 69, ('T', 'W'): 128, ('T', 'Y'): 92,
('V', 'A'): 64, ('V', 'C'): 192, ('V', 'D'): 152, ('V', 'E'): 121, ('V', 'F'): 50, ('V', 'G'): 109, ('V', 'H'): 84, ('V', 'I'): 29, ('V', 'K'): 97, ('V', 'L'): 32, ('V', 'M'): 21, ('V', 'N'): 133, ('V', 'P'): 68, ('V', 'Q'): 96, ('V', 'R'): 96, ('V', 'S'): 124, ('V', 'T'): 69, ('V', 'V'): 0, ('V', 'W'): 88, ('V', 'Y'): 55,
('W', 'A'): 148, ('W', 'C'): 215, ('W', 'D'): 181, ('W', 'E'): 152, ('W', 'F'): 40, ('W', 'G'): 184, ('W', 'H'): 115, ('W', 'I'): 61, ('W', 'K'): 110, ('W', 'L'): 61, ('W', 'M'): 67, ('W', 'N'): 174, ('W', 'P'): 147, ('W', 'Q'): 130, ('W', 'R'): 101, ('W', 'S'): 177, ('W', 'T'): 128, ('W', 'V'): 88, ('W', 'W'): 0, ('W', 'Y'): 37,
('Y', 'A'): 112, ('Y', 'C'): 194, ('Y', 'D'): 160, ('Y', 'E'): 122, ('Y', 'F'): 22, ('Y', 'G'): 147, ('Y', 'H'): 83, ('Y', 'I'): 33, ('Y', 'K'): 85, ('Y', 'L'): 36, ('Y', 'M'): 36, ('Y', 'N'): 143, ('Y', 'P'): 110, ('Y', 'Q'): 99, ('Y', 'R'): 77, ('Y', 'S'): 144, ('Y', 'T'): 92, ('Y', 'V'): 55, ('Y', 'W'): 37, ('Y', 'Y'): 0}


def get_aa_single(long_form):
    '''Convert from long form amino acid names (HIS, Histidine, etc) to single letter IDs (H) either strings or lists of strings'''
    
    if type(long_form)==type([]):  #if list	
        return([aa_dict[i.replace(' ', '').strip().upper()] for i in long_form])
    else:
        return aa_dict[long_form.replace(' ', '').strip().upper()]

def get_aa_mut(dnaseq, pos, alt_allele):
    '''	dnaseq: CDS sequence, codon aligned
        pos: 1-indexed position of mutated locus
        alt_allele: alternate allele at the position
        
        returns: WT AA, mutant AA, missense/nonsense/synonymous label '''
        
    pos -= 1
    seq = list(dnaseq.upper())
    alt_allele = alt_allele.upper()
    
    codon_pos = pos%3
    wt_codon = seq[pos-codon_pos: pos-codon_pos+3]
    mut_codon = wt_codon[:]
    mut_codon[codon_pos] = alt_allele
    
    wt_aa = codon_dict[''.join(wt_codon)]
    mut_aa = codon_dict[''.join(mut_codon)]
    if wt_aa == mut_aa:
        label = 'synonymous'
    elif mut_aa == '*':
        label = 'nonsense'
    else:
        label = 'missense'
        
    return wt_aa, mut_aa, label

###############################################################################
################################### UNIPROT ###################################
###############################################################################

def fetch_uniprot(uniprots, attributes):
    '''Fetch information about a UniProt ID from uniprot.org
    UniProt attributes can be a single string, which will return a single string value,
    or a list of attributes, which will return a dictionary of results attribute => value
    
    accepted attributes:
    
    citation | clusters | comments | domains | domain | ec | id
    entry name | existence | families | features | genes | go
    go-id | interactor | keywords | last-modified | length
    organism | organism-id | pathway | protein names | reviewed
    sequence | 3d | version | virus hosts
    
    For more info, see: http://www.uniprot.org/help/programmatic_access
    '''
    
    if type(uniprots) == str:
        uniprots = [uniprots,]
    
    if type(attributes) == str:
        attributes = [attributes,]
    
    if 'id' not in attributes:
        attributes.append('id')
    
    # EDIT: Shayne Wierbowski July 24, 2018
    # UniProt API seems to have changed, joining uniprot IDs by commas no longer seems to work
    # I have modified this code to join uniprot IDs by "+OR+" instead as if you were building a query
    # directly as described at https://www.uniprot.org/help/api_queries
    # I am not sure why a similar change was not required to the columns field
    # Changed from...
    #uniprots = ','.join(uniprots)
    # To...
    uniprots = '+OR+'.join(uniprots)
    columns = ','.join(attributes)
    
    #build and send API request
    # EDIT: Shayne Wierbowski July 24, 2018
    # UniProt API now seems to require https rather than http url
    #Changed from...
    #url = 'http://www.uniprot.org/uniprot/' by Shayne Wierbowski July 24, 2018 because
    # To...
    url = 'https://www.uniprot.org/uniprot/'
    params = {'query': uniprots,
              'format': 'tab',
              'columns': columns}
    data = urllib.parse.urlencode(params)
    request = urllib.request.Request(url, data)
    request.add_header('User-Agent', 'Python mjm659@cornell.edu')
    response = urlopen(request)
    data = response.readlines()
    ###

    if data==[] or len(data)<2:
        return {}
    
    # EDITS MADE BY SHAYNE WIERBOWSKI ON OCTOBER 30m 2019
    # Added extra_dict / warnings for extra IDs that get fetched but are not
    # requested. Also added handling for "Merged Into" rows (e.g.)
    #
    # P18463	1B37_HUMAN							Merged into P01889.	
    #
    return_dict = {}
    extra_dict = {}
    redo = []
    
    for line in data[1:]:
        
        line_data = line.strip().split('\t')
        if len(line_data) != len(attributes):
            if("Merged into" in line):
                redo.append((line.strip().split("\t")[0], line.strip().split("\t")[-1].split()[-1].replace(".", "")))
            continue
        
        entry_dict = {}
        for i, a in enumerate(attributes):
            entry_dict[a] = line_data[i]
        if entry_dict['id'] in uniprots:
            return_dict[entry_dict['id']] = entry_dict
        else:
            extra_dict[entry_dict['id']] = entry_dict
    
    if(len(extra_dict) != 0):
        # Apparrently this just happens all the time, but I can't figure out why?
        #warnings.warn("Unrequested UniProts ({0}) were returned. This should not happen.".format(";".join(extra_dict.keys())))
        pass
    
    for p1, p2 in redo:
        if(p2 in return_dict):
            return_dict[p1] = copy.deepcopy(return_dict[p2])
            return_dict[p1]["id"] = p1
        elif(p2 in extra_dict):
            return_dict[p1] = copy.deepcopy(extra_dict[p2])
            return_dict[p1]["id"] = p1
        warnings.warn("Requested UniProt ID, {0}, appears to have been replaced by {1}".format(p1, p2))
    
    # END EDITS
    
    return return_dict
    
    
def fetch_uniprot_sequence(uniprot, uniprotDict={}):
    '''Fetch sequence for a UniProt ID from uniprot.org
    Optional - Pass a dictionary mapping uniprot IDs to sequence. Will only call the web service if uniprot is not in dict.'''
    
    if uniprot in uniprotDict:
        #retrieve uniprot from provided dictionary
        return uniprotDict[uniprot]
    else:
        #retrieve uniprot from the web
        try:
            return fetch_uniprot(uniprot, 'sequence')[uniprot]['sequence']
        except KeyError:
            return ''
    

def convert_to_uniprot(query):
    '''Return UniProt results for any identifier.'''
    
    query = query.replace(' ', '%20')
    
    columns = 'id,reviewed,organism,organism-id,genes,protein%20names,length,sequence'
    link = 'http://www.uniprot.org/uniprot/?query=%s&sort=score&format=tab&columns=%s' %(query, columns)
    data = urlopen(link).read()
    
    if data == '':
        return None
    
    header = data.split('\n')[0].split('\t')
    results = [r.split('\t') for r in data.split('\n')[1:-1]]

    return_list = []
    
    for r in results:
        cur_dict = {}
        for i in range(len(header)):
            cur_dict[header[i]] = r[i]
        return_list.append(cur_dict)
        
    return return_list


def parse_uniprot(uniprot_acc):
    '''fetch common data fields from a UniProt file.
    
    - uniprot_acc can be either a local file (if path exists), or a UniProt ID which will fetch data from the web.
    
    Parse common fields out of a UniProt file.'''
    
    if os.path.exists(uniprot_acc):
        file_handle = open(uniprot_acc)
    else:
        file_handle = urlopen('http://www.uniprot.org/uniprot/%s.txt' %(uniprot_acc.upper()))
        

    data =   {'synonyms': [],
            'reviewed': False,
            'accession': '',
            'taxid': '',
            'length': '',
            'aaseq': '',
            'recName': '',
            'entrez': set(),
            'ucsc': [],
            'geneName': [],
            'dbSNP': [],
            'refT': [],
            'refP': [],
            'ensT': [],
            'ensP': [],
            'ensG': [],
            }
    
    #~ data['uniprot'] = path.splitext(path.basename(f))[0]   #CHANGE THIS
    seq_start = False
    
    for l in file_handle:
        
        if seq_start == True:
            data['aaseq'] += l.strip().replace(' ', '').replace('//', '')
        
        if l[:5]=='AC   ':
            data['synonyms'] = [up.strip() for up in l[5:].strip().split(';') if len(up.strip())>1]
        
        if l[:5]=='ID   ':
            if 'Reviewed;' in l:
                data['reviewed'] = True
        
        if l[:5]=='ID   ':
            data['accession'] = l[5:].strip().split()[0]
        
        if l[:16]=='OX   NCBI_TaxID=':
            data['taxid'] = l[16:].split(';')[0]
        
        if l[:13]=='DR   GeneID; ':
            rs = [i.replace(';', '') for i in l.replace('DR   GeneID; ', '').strip()[:-1].split()]
            for r in rs:
                if r != '-': data['entrez'].add(r)
        
        if l[:13]=='DR   RefSeq; ':
            line = l.replace('DR   RefSeq; ', '').strip()
            if len(findall('\[.*?\-.*?\]', line)) > 0 and '-1]' not in line:   #don't keep non-canonical isoform matches
                continue
            ids = findall('N[M,P]_\d*\.\d', line)
            #~ print ids
            for e in ids:
                if e[:3]=='NM_': data['refT'].append(e)
                if e[:3]=='NP_': data['refP'].append(e)
                
                
        if l[:11]=='DR   UCSC; ':
            line = l.replace('DR   UCSC; ', '').strip()
            ucsc_id, organism = [r.strip() for r in line.split(';')]
            if len(findall('\[.*?\-.*?\]', organism)) > 0 and '-1]' not in organism:   #don't keep non-canonical isoform matches
                continue
            data['ucsc'].append(ucsc_id)
                
        if l[:14]=='DR   Ensembl; ':
            line = l.replace('DR   Ensembl; ', '').strip()
            if len(findall('\[.*?\-.*?\]', line)) > 0 and '-1]' not in line:   #don't keep non-canonical isoform matches
                continue
            ids = findall('ENS[T,G,P]\d*', line)
            for e in ids:
                if e[:4]=='ENST': data['ensT'].append(e)
                if e[:4]=='ENSP': data['ensP'].append(e)
                if e[:4]=='ENSG': data['ensG'].append(e)
        
        if l[:10]=='GN   Name=':
            rs = [i.strip() for i in l.replace('GN   Name=', '').replace('Synonyms=', '').replace('ORFNames=', '').replace('/', ',').replace(';', ',').strip()[:-1].split(',')]
            for r in rs:
                if r != '-': data['geneName'].append(sub('\{.*?\}', '', r).strip())
                
        if l[:19]=='DE   RecName: Full=':
            data['recName'] = l.replace('DE   RecName: Full=', '').strip()[:-1]
        
        if l[:12]=='FT   VARIANT' and 'dbSNP' in l:
            tokens = l.strip().split()
            if tokens[5] != '->': continue
            if len(tokens) != 9: continue
            
            
            pos1 = tokens[2]; pos2 = tokens[3]
            if pos1 != pos2:
                print(uniprot, l)
                exit()
            
            orig_aa = tokens[4]
            new_aa = tokens[6]
            dbsnp_ref = tokens[-1].split(':')[-1].replace(').', '')
            
            data['dbSNP'].append('%s%s%s(%s)' %(orig_aa, pos1, new_aa, dbsnp_ref))
            
        
        if l[:13]=='SQ   SEQUENCE':
            data['length'] = l.replace('SQ   SEQUENCE', '').strip().split(' AA;')[0]
            seq_start = True
                    
    file_handle.close()
    
    return data


###############################################################################
###################################   PDB   ###################################
###############################################################################

def is_binary_file(filename):
    ''' check if file is binary (adapted from http://stackoverflow.com/questions/898669/how-can-i-detect-if-a-file-is-binary-non-text-in-python) '''
    textchars = bytearray([7,8,9,10,12,13,27]) + bytearray(range(0x20, 0x100))
    return bool(open(filename, 'rb').read(1024).translate(None, textchars))


def natural_keys(text):
    '''natural sorting (adapted from http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside)
          Useful for sorting PDB residues that have letters in them
       
       example:
       
       >> from mjm_tools import natural_keys
    
       >> my_list = ['1', '3', '2', '7', '2B']   
       >> my_list.sort(key=natural_keys)
       ['1', '2', '2B', '3', '7']
          
          '''
    def atoi(text): return int(text) if text.isdigit() else text
    return [ atoi(c) for c in regexsplit('(\d+)', text) ]


def dist3D(point1, point2):
    '''Find the 3D distance beteen two points as tuples of floats. This implementation is much faster than numpy's linalg.norm'''
    return ((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2 + (point1[2]-point2[2])**2)**0.5	


def open_pdb(structure, verbose=True, try_web=True):
    '''Return an opened PDB file handle from STDIN, file, local PDB cache, or web'''

    # STDIN
    if "<open file '<stdin>', mode 'r' at" in str(structure):
        pdb_filehandle = structure

    # AS UNCOMPRESSED PDB FILE
    elif os.path.exists(structure) and is_binary_file(structure)==False:   #file exists and is a text-based file
        pdb_filehandle = open(structure, 'r')
        
    # AS GZIPPED PDB FILE	
    elif os.path.exists(structure) and is_binary_file(structure)==True:    #file exists and is likely a gzipped file
        try:
            testopen = gzopen(structure, 'r')
            testopen.readline()
            testopen.close()
            pdb_filehandle = gzopen(structure, 'r')
        except IOError:
            if(verbose):
                print('Invalid structure file-type. Structure file must be a plain-text PDB file or a gzipped PDB file.')
            return

    # AS PDB FILE FROM LOCAL COPY OF THE PDB -OR- FROM THE WEB
    elif len(structure)==4:
        
        pdb_storage_path = os.path.join(PDB_DATA_DIR, '%s/pdb%s.ent.gz' %(structure[1:3].lower(), structure.lower()))
        
        #local file
        if os.path.exists(pdb_storage_path):
            pdb_filehandle = gzopen(pdb_storage_path, 'r')
        #try the web
        elif(try_web):
            try:
                pdb_filehandle = urlopen('http://www.rcsb.org/pdb/files/%s.pdb' %(structure.upper()))
            except HTTPError:
                if(verbose):
                    print('Invalid structure input: %s. Not found as local file, as PDB structure in %s, or on the web.' %(structure, PDB_DATA_DIR))
                return	
        else:
            return
    else:
        if(verbose):
            print('Invalid structure input: %s. Not found as local file, and wrong number of characters for direct PDB reference.' %(structure))
        return
    
    return pdb_filehandle


def rotate_pdb(structure, out_file, radians=1.3):
    '''Returns a rotated version of the original PDB file, with all coordinates updated to new rotated positions
    while maintaining inter-atom spatial relationships. Rotates about the center of mass by specified number of radians.
    1.3 radians by default ~ 75 degrees'''
    
    def rotation_matrix(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians. FROM: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
        """
        axis = np.asarray(axis)
        theta = np.asarray(theta)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2)
        b, c, d = -axis*math.sin(theta/2)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
    
    #-----------------------------------------------------------------------------

    pdbinfile = open_pdb(structure)
    
    out_lines = []
    coordinates = []
    atomline2coords = {}      # full pdb line: {'original': (x,y,z), 'transformed': (x,y,z)}
    
    recordTypes = set(['ATOM', 'HETATM'])
    
    for l in pdbinfile:
        recordName = l[:6].strip()
        
        out_lines.append(l)
        
        if len(l) < 54 or recordName not in recordTypes:
            continue
        
        x = float(l[30:38].strip())
        y = float(l[38:46].strip())
        z = float(l[46:54].strip())
        
        atomline2coords[l] = (x,y,z)
        
    #calculate center of structure:
    x_mean = np.mean([atomline2coords[l][0] for l in atomline2coords])
    y_mean = np.mean([atomline2coords[l][1] for l in atomline2coords])
    z_mean = np.mean([atomline2coords[l][2] for l in atomline2coords])
    
    rotation_axis = (x_mean, y_mean, z_mean)
    theta = 1.2
    rot_mat = rotation_matrix(rotation_axis,theta)
    
    output = open(out_file, 'w')
    
    for l in out_lines:
        if l not in atomline2coords:
            output.write(l)
            continue
        
        new_coords = np.dot(rot_mat, atomline2coords[l])
        
        output.write(l[:30] + str(round(new_coords[0],3)).rjust(8) + str(round(new_coords[1],3)).rjust(8) + str(round(new_coords[2],3)).rjust(8) + l[54:])
    
    output.close()


def parse_naccess(naccess_output):
    '''	Parse naccess output into dictionary.
    
        Example head of naccess output file (.rsa):
    
        REM  Relative accessibilites read from external file "/usr/local/naccess2.1.1/standard.data"
        REM  File of summed (Sum) and % (per.) accessibilities for 
        REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
        REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
        RES SER A  14    28.88  24.8  22.22  28.5   6.66  17.3  18.69  38.5  10.19  15.0
        RES ALA A  15    10.35   9.6   0.00   0.0  10.35  26.8   2.33   3.3   8.01  21.9
        RES ASN A  16    37.02  25.7  36.33  34.2   0.69   1.8  25.49  55.1  11.53  11.8
        RES LEU A  17   130.51  73.1 105.12  74.5  25.39  67.7 105.51  74.1  25.00  68.8
        RES ASP A  18   103.34  73.6 100.05  97.4   3.29   8.7  33.55  68.1  69.79  76.6
        RES HIS A  19   143.10  78.2 123.76  84.1  19.34  54.0  75.48  77.7  67.62  78.9
        
    '''
    
    parsed_output = []
    
    for l in naccess_output:
        if l[:3] != 'RES':
            continue
        
        res, chain, res_num = l[3:7].strip(), l[7:9].strip(), l[9:14].strip()
        combined_key = l[4:14]
        
        all_atoms_abs, all_atoms_rel = float(l[14:22].strip()), float(l[22:28].strip())
        total_side_abs, total_side_rel = float(l[28:35].strip()), float(l[35:41].strip())
        main_chain_abs, main_chain_rel = float(l[41:48].strip()), float(l[48:54].strip())
        non_polar_abs, non_polar_rel = float(l[54:61].strip()), float(l[61:67].strip())
        all_polar_abs, all_polar_rel = float(l[67:74].strip()), float(l[74:80].strip())
        
        parsed_output.append({	'res': res, 'chain': chain, 'res_num': res_num, 'combined_key': combined_key,
                            'all_atoms_abs': all_atoms_abs, 'all_atoms_rel': all_atoms_rel,
                            'total_side_abs': total_side_abs, 'total_side_rel': total_side_rel,
                            'main_chain_abs': main_chain_abs, 'main_chain_rel': main_chain_rel,
                            'non_polar_abs': non_polar_abs, 'non_polar_rel': non_polar_rel,
                            'all_polar_abs': all_polar_abs, 'all_polar_rel': all_polar_rel })
        
    return parsed_output
        
    
    
def naccess(pdb_file, parsed=False):
    # Added 2021_07_08 by Shayne Wierbowski
    # Naccess seems to break when provided with relative paths
    pdb_file = os.path.realpath(pdb_file)
    '''Run NACCESS and return the results.'''
    cwd = os.getcwd()     #naccess writes to the current working directory, so save current directory, and move to scatchDir to write output files
    os.chdir(os.path.dirname(pdb_file))   #write naccess output files to directory of PDB file

    raw_naccess_output = []
    
    # Delete output files if they already exist (EDIT made by Shayne July 12, 2018)
    # Leaving these files can be a problem when running naccess multiple times (i.e. irescalc.py)
    # Since if the naccess command fails, this function will see the previous outputs and not detect the failure
    if(os.path.exists(os.path.splitext(pdb_file)[0]+'.rsa')):
        os.remove(os.path.splitext(pdb_file)[0]+'.rsa')
    if(os.path.exists(os.path.splitext(pdb_file)[0]+'.asa')):
        os.remove(os.path.splitext(pdb_file)[0]+'.asa')
    
    _, _ = Popen([NACCESS_PATH, pdb_file], stdout=PIPE, stderr=PIPE).communicate()
    try:
        raw_naccess_output += open(os.path.splitext(pdb_file)[0]+'.rsa').readlines()
        
        # Make sure output files are not empty (EDIT made by Shayne July 12, 2018)
        # See above
        if(os.stat(os.path.splitext(pdb_file)[0]+'.rsa').st_size == 0):
            raise IOError('ERROR: Naccess .rsa file is empty. We suspect this is an edge case where Naccess cannot calculate ASA for extremely large chains.'
                          ' The following command was attempted: %s %s' %(NACCESS_PATH, pdb_file))
            # print('WARNING: This error message is only a temporary measure. The longterm solution will involve identifying an alternative SASA calculator. If you see this message, please contact Shayne.')
            # exit()
    except IOError:
        raise IOError('ERROR: Naccess .rsa file was not written. The following command was attempted: %s %s' %(NACCESS_PATH, pdb_file))
        exit()

    os.chdir(cwd)  #move back to old cwd
    
    if parsed:
        return parse_naccess(raw_naccess_output)
    else:
        return raw_naccess_output



def freesasa(pdb_file, parsed=False, warn=True):
    if(warn):
        print ("Hello. You've found a secret warning message left by Shayne. "
               "FreeSASA cannot calculate relative side-chain accessibilities for GLY (has no side chain) "
               "and instead reports 'N/A' in its output. By contrast, NACCESS (which we used to use) include "
               "the alpha-carbon in the RSA calculation (so this was not a problem). For our current purposes, "
               "we never use these values (we only use the all-atom accessibility columns). Having the 'N/A' in"
               " the output breaks Michael's parser for RSA output files, so I modify these values to '0.0' "
               "instead. If you happen to be using this function and are interested in the raw relative side-chain "
               "accessibilities be aware of this problem; there should be no valid value for GLY, but they are all "
               "reported as 0.0 by this function). As long as you understand the implications of this message and "
               "it does not affect your analysis, set the 'warn' parameter of this funciton to 'False' to disable this message.")
    
    '''Run FREESASA and return the results.'''
    cwd = os.getcwd()     #naccess writes to the current working directory, so save current directory, and move to scatchDir to write output files
    os.chdir(os.path.dirname(pdb_file))   #write naccess output files to directory of PDB file
    
    raw_naccess_output = []
    
    # Delete output files if they already exist (EDIT made by Shayne July 12, 2018)
    # Leaving these files can be a problem when running naccess multiple times (i.e. irescalc.py)
    # Since if the naccess command fails, this function will see the previous outputs and not detect the failure
    if(os.path.exists(os.path.splitext(pdb_file)[0]+'.rsa')):
        os.remove(os.path.splitext(pdb_file)[0]+'.rsa')
    if(os.path.exists(os.path.splitext(pdb_file)[0]+'.asa')):
        os.remove(os.path.splitext(pdb_file)[0]+'.asa')
    
    _, _ = Popen([FREESASA_PATH, "--format=rsa", pdb_file], stdout=open(os.path.splitext(pdb_file)[0]+'.rsa', "w+"), stderr=PIPE).communicate()
    try:
        # We have to manually replace the N/A that appear for the side-chain rel SASA in Glycine (has no side chain) with 0.0 as a placeholder
        # NOTE: This is technically an error / shortcoming with FreeSASA comapred with NACCESS, but for our purposes we don't actually use the
        # side-chain rel SASA columns, only the all-atoms rel SASA columns. This COULD be an issue if someone else ever uses raw SASA values
        # from FreeSASA for their purposes
        raw_naccess_output += [l.replace("N/A", "0.0") for l in open(os.path.splitext(pdb_file)[0]+'.rsa').readlines()]
        
        # Make sure output files are not empty (EDIT made by Shayne July 12, 2018)
        # See above
        if(os.stat(os.path.splitext(pdb_file)[0]+'.rsa').st_size == 0):
            print('ERROR: FreeSASA .rsa file is empty. I have never encountered this error before. May still be related to using extremely large chains. The following command was attempted: %s --format=rsa %s' %(FREESASA_PATH, pdb_file))
            print('WARNING: This error message is only a temporary measure. The longterm solution will involve identifying an alternative SASA calculator. If you see this message, please contact Shayne.')
            exit()
    except IOError:
        raise IOError('ERROR: FreeSASA .rsa file was not written. The following command was attempted: %s --format=rsa %s' %(FREESASA_PATH, pdb_file))
        exit()
    
    os.chdir(cwd)  #move back to old cwd
    
    if parsed:
        return parse_naccess(raw_naccess_output)
    else:
        return raw_naccess_output



def get_chains(structure):
    '''Returns a list of chains in a PDB file'''
    
    pdbinfile = open_pdb(structure)
    
    chains = []
    chains_set = set()  #also build set alongside list for fast lookup of whether we've encountered a chain yet
    record_types=set(['ATOM', 'HETATM'])
    
    for l in pdbinfile:
        recordName = l[:6].strip()
        
        if recordName not in record_types:
            continue
            
        chainID = l[21]
        
        if chainID not in chains_set:
            chains.append(chainID)
            chains_set.add(chainID)
        
    return chains
    

def get_residues(structure, chain):
    '''Returns a list of residues in a specified PDB file and chain'''
    
    pdbinfile = open_pdb(structure)
    
    residues = []
    residues_set = set()  #also build set alongside list for fast lookup of whether we've encountered a residue yet
    
    for l in pdbinfile:
        
        if len(l) < 30: continue
        
        chainID = l[21]
        resSeq = l[22:27].strip()
        
        if chainID != chain:
            continue
        
        if resSeq not in residues_set:
            residues.append(resSeq)
            residues_set.add(resSeq)
        
    return residues_set


def get_chainseqs(structure, chain=None):
    '''Returns a string sequence of residues in a given structure and chain'''
    
    pdbinfile = open_pdb(structure)
    
    chain2sequence = defaultdict(str)
    chain2residues_set = defaultdict(set)  #also build set alongside list for fast lookup of whether we've encountered a residue yet
    
    record_types = set(['ATOM', 'HETATM'])
    
    for l in pdbinfile:
        
        recordName = l[:6].strip()
        
        if recordName not in record_types:
            continue
        
        if len(l) < 30: continue
        
        resName = l[17:20].strip()
        chainID = l[21]
        resSeq = l[22:27].strip()
        
        if chain!= None and chainID != chain:
            continue
        
        if resSeq not in chain2residues_set[chainID]:
            chain2sequence[chainID] += get_aa_single(resName)
            chain2residues_set[chainID].add(resSeq)
    
    if chain!=None:
        return chain2sequence[chain]
    else:
        return dict(chain2sequence)

def profit_overlap(res1, res2, chain, profit_prefix):
    '''Return a profit-formatted zone struct of overlapping residues
    profit_prefix = zone, rzone, etc.'''
        
    overlap_res = res1.intersection(res2)
    
    all_res = zip_res_range(sorted(list(overlap_res), key=natural_keys))
    
    profit_zones = ''
    
    for r in all_res[1:-1].split(','):
        profit_zones += profit_prefix + ' ' + chain + r.split('-')[0] + '-' + chain + r.split('-')[-1] + '\n'
    
    return profit_zones.strip(), len(overlap_res)
    


def xtool(tool, *args, **kwargs):
    ''' Python function front-end to srescalc.py, irescalc.py, crescalc.py, etc from MJM_TOOLS folder
    *args should contain positional arguments, **kwargs contains flagged option arguments'''
    
    tool_path = os.path.join(TOOLS_PATH, tool)
    command = [tool_path,] + list(args)
    for k, v in kwargs.items(): command += ['-%s' %(k), str(v)]
    out, err = Popen(command, stderr=PIPE, stdout=PIPE).communicate()
    return out


# EDIT MADE BY SHAYNE WIERBOWSKI 2019_06_24
# - Added the only_first_model flag so that the atom extraction finishes after
#   the first model and does not include duplicate atoms
# - The way the function was set up, the PDB metadata to separate different
#   models was lost and so this function sometimes created broken PDBs.
def extract_atoms(structure, out_file, chain_dict = {}, chain_map = {}, atoms = set(), record_types=set(['ATOM', 'HETATM']), write_end=False, fix_modbase=False, only_first_model=False):
    ''' Extract only specific residues from specific chains
    and specific atoms (i.e. backbone) within those residues.
    
    -structure can be file, gzipped_file, or pdbID (read from resources)
    -out_file is location to write parsed .pdb output file
    
    -chain_dict is dictionary containing chains and residues to take from input file
        i.e. chain_dict = { 'A': set(['1','2','3','4','5','10','11','11B','12'...]), 'B': set([...]) }
     OR chain_dict can contain mapppings to new residue names: (i.e for mapping to UniProt residues)
        i.e. chain_dict = { 'A': {'1':'27', '2':'28', '3':'29', '4':'30', '5':'31'},  'B': {...} }
        
    -chain_map is a dictionary mapping native chains to new chain IDs in output_file
        i.e. chain_map = {'G': 'A', 'C': 'B'}  maps all chain G to chain A and all chain C to chain B, resulting chains will be alphabetically stored in output file
    
    -atoms is a set containing atom types to keep in final structure
        i.e. atoms = set(['CA', 'C', 'N', O, ...])
    
    -record_types are the types of records kept (field 1 - i.e. ATOM, HETATM, ANISOU, etc.)
    
    Empty vars mean take all--i.e. chain_dict = {} means take all chains
                                   chain_dict = {'A': set(), 'B': set()} means take all residues in chains A and B
                                   atoms = set() means take all atoms in whichever residues and chains indicated in chain_dict
                                   record_types = set() means take all record types
                                   
    -write_end = True writes 'END' as the last line of the new file            
    
    -fix_modbase removes problematic fields that often show up in modbase models.                    
                                   '''
    
    #-----Open input structure-----
    
    pdbinfile = open_pdb(structure)
    
    # EDIT Shayne Wierbowski 2020_11_27
    # - Added check and error reporting if open_pdb returns None
    if(pdbinfile is None):
        raise TypeError("Invalid Structure ({0})".format(structure))
    
    #-----Parse input structure-----
    
    new_lines = defaultdict(list)
    integers = set(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])
    
    for l in pdbinfile:
        #pdb fields described here: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
        recordName = l[:6].strip()

        if record_types != set() and recordName not in record_types:
            # EDIT MADE BY SHAYNE WIERBOWSKI 2019_06_24
            # Stop after the end of the first model if optional flag is set (False by default)
            if(only_first_model and recordName == "ENDMDL"):
                break
            # END EDIT
            continue
            
        atomName = l[12:16].strip()
        chainID = l[21]
        resName = l[17:20].strip()
        resSeq = l[22:27].strip()     # resSeq 22:26, and 26:27 is the insertion character for special residues, i.e. 37A
            
        new_line = l
        
        if chainID not in chain_dict and chain_dict != {}:  #if chain not found within a set of specified chains. (chain could not be found, but the user may not have indicated ANY chains, in which case we use ALL chains, not none.)
            continue
        
        if chain_dict != {} and len(chain_dict[chainID]) != 0:  #if list of residues is NOT an empty structure (meaning take all residues), check that current residue is one we want to save
            if resSeq not in chain_dict[chainID]:
                continue
            if type(chain_dict[chainID])==dict: #mapping to new residue names
                
                if chain_dict[chainID][resSeq][-1] in integers:
                    new_resSeq = chain_dict[chainID][resSeq].rjust(4) + ' '
                else:
                    new_resSeq = chain_dict[chainID][resSeq].rjust(5)
                    
                new_line = new_line[:22] + new_resSeq + new_line[27:]
                
        if atoms != set() and atomName not in atoms:
            continue
            
        if fix_modbase:
            new_line = new_line[:73] + ' '*4 + new_line[13] + ' '*3 + '\n'  
        
        if chainID in chain_map:
            new_line = new_line[:21] + chain_map[chainID] + new_line[22:]
            index_chain = chain_map[chainID]
        else:
            index_chain = chainID
        
        new_lines[index_chain].append(new_line)
    
    #-----Write output structure file-----
        
    output = open(out_file, 'w')
    for chain in sorted(new_lines.keys()):
        for line in new_lines[chain]:
            output.write(line)
    if write_end: output.write('END\n')
    output.close()


# EDIT MADE BY SHAYNE WIERBOWSKI 2019_02_20
#
# Added a saved global sifts_data file so that when running this
# function repeatedly it does not need to reparse the entire
# sifts file every time.
#
# Also added a parameter that can be turend off to ignore this
# behavior.
sifts_cache = [None]

# EDIT MADE BY SHAYNE WIERBOWSKI 2019_06_24
# - Added only_first_model parameter
def pdb_sifts_map(pdb, chain, out_file, cache_sifts_data=True, only_first_model=False):
    '''Take in pdb ID and Chain, output file with only UniProt residues mapped to UniProt indices'''
    
    pdb = pdb.upper()
    
    # EDITED BY SHAYNE WIERBOWSKI 2019_02_20 to improve efficiecy when being called many times
    if(cache_sifts_data and sifts_cache[0] == None):
        sifts_data = parse_dictionary_list(SIFTS_MAPPING_FILE)
        sifts_cache[0] = sifts_data
    elif(cache_sifts_data):
        sifts_data = sifts_cache[0]
    else:
        sifts_data = parse_dictionary_list(SIFTS_MAPPING_FILE)		
    # END EDITS
    
    for d in sifts_data:
        if d['PDB'] == pdb and d['Chain']==chain:
            uniprot_res = unzip_res_range(d['MappableResInPDBChainOnUniprotBasis'])
            pdb_res = unzip_res_range(d['MappableResInPDBChainOnPDBBasis'])
            break
    else:
        print('No SIFTS entry found for given PDB')
        return
    
    residue_mapping_dict = dict([(pdb_res[i], uniprot_res[i]) for i in range(len(uniprot_res))])
    residue_mapping_dict = {chain: residue_mapping_dict}
    
    # EDIT MADE BY SHAYNE WIERBOWSKI 2019_06_24
    # Added only_first_model parameter
    extract_atoms(pdb, out_file, chain_dict = residue_mapping_dict, only_first_model=only_first_model) 
    
    
    


def pdb_segid_to_chain(infile, outfile):
    '''Replaces haddock script tools/pdb_segid-to-chain which copies the segid given in column 73 to the chainID in column 22
        Also, replace BB with CA to account for sed command in original haddock i-rmsd script'''

    output = open(outfile, 'w')
    for l in open(infile):
        if l[:4]=='ATOM' or l[:6]=='HETATM' and l[72] != ' ':
            new_line = l[:21] + l[72] + l[22:]
        else:
            new_line = l
        output.write(new_line.replace('BB', 'CA'))
    output.close()


###############################################################################
###################################  OMIM  ####################################
###############################################################################

def to_ascii(any_str):
    '''convert any ascii or unicode string to ascii'''
    return unicodedata.normalize('NFKD', unicode(any_str, "utf-8")).encode('ascii', 'ignore')
    

def fetch_omimID(search_term):
    '''Finds top search hit on omim.org for given search term. Returns disease name and OMIM ID.'''
    
    omimID, omim_disease = 'n/a', 'n/a'
    
    try:
        search_term = to_ascii(search_term)
        url = 'http://omim.org/search?index=entry&search='+search_term.replace(' ', '+')+'&sort=score+desc%2C+prefix_sort+desc&start=1&limit=20&format=tab'
        request = urllib.request.Request(url, headers={'User-Agent' : 'Mozilla/5.0'}) 
        data = [l.strip() for l in urllib.request.urlopen(request).readlines()]
            
        for e in data:
            if len(e) > 1 and (e[0] == '#' or e[0] == '%'):
                omimID = e[1:7]
                omim_disease = e.split('\t')[1]
                break
    except:
        pass
            
    return omimID, omim_disease

###############################################################################
###################################  OTHER  ###################################
###############################################################################

def time_elapsed(start_time):
    elapsed = time.time()-start_time

    if elapsed < 60: # minute
        return '%.2f sec' %(elapsed)
    elif elapsed < 3600: #hour
        return '%i min %i sec' %((elapsed%3600)//60, (elapsed%3600)%60)
    elif elapsed < 86400: #day
        return '%i hr %i min %i sec' %(elapsed//3600, (elapsed%3600)//60, (elapsed%3600)%60)
    elif elapsed < 172800: #2 days
        return '%i day %i hr %i min %i sec' %(elapsed//86400, (elapsed%86400)//3600, (elapsed%3600)//60, (elapsed%3600)%60)
    else:
        return '%i days %i hr %i min %i sec' %(elapsed//86400, (elapsed%86400)//3600, (elapsed%3600)//60, (elapsed%3600)%60)

import subprocess, threading

class TimeoutCmd(threading.Thread):
    '''Allows running a shell command with a timeout threshold, afterwhich the command will stop execution even if it has not finished.
    
    Example usage:
      from mjm_tools import TimeoutCmd
      timeout = 3
      myCmd = TimeoutCmd(['blastp', '-query', 'Legionella_selected.faa', '-db', 'Legionella_all.faa', '-outfmt', '6'], timeout)
      myCmd.Run()
      print 'output', myCmd.out
    '''
    
    def __init__(self, cmd, timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout
    
    def run(self):
        self.p = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE)
        self.out, self.err = self.p.communicate()
        
    def Run(self):
        self.start()
        self.join(self.timeout)
        
        if self.is_alive():
            self.p.terminate()      #use self.p.kill() if process needs a kill -9
            self.join()
            self.out = 'FAILED TO FINISH'
            


# import smtplib
# from email.MIMEMultipart import MIMEMultipart
# from email.MIMEBase import MIMEBase
# from email.MIMEText import MIMEText
# from email import Encoders

# def mail(to, subject, text):
#     '''send an email with the hyulab gmail account'''
    
#     gmail_usr = "hyulab@gmail.com"
#     gmail_pwd = "W0"+"drkj"[1:-1]+"h;0;r;d!".replace(';', '')
    
#     msg = MIMEMultipart()

#     msg['From'] = gmail_usr
#     msg['To'] = to
#     msg['Subject'] = subject
    
#     msg.attach(MIMEText(text))

#     mailServer = smtplib.SMTP("smtp.gmail.com", 587)
#     mailServer.ehlo()
#     mailServer.starttls()
#     mailServer.ehlo()
#     mailServer.login(gmail_usr, gmail_pwd)
#     mailServer.sendmail(gmail_usr, to, msg.as_string())
#     # Should be mailServer.quit(), but that crashes...
#     mailServer.close()
