"""
Determine triple/double correlated sources               
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu      

Use the cross correlate variable sources across the 3 EPIC cameras
from the previous step and organise the results into triple matches
and double matches. Sources which are not cross correlated are
preseved and labled as singles.
"""
import argparse
import re

import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.table import Table

from exodus_utils import check_correlation


###
# Parsing arguments
###
parser = argparse.ArgumentParser()

# Variability parameters
parser.add_argument("-path", help="Path to the folder containing the observation files", type=str)

parser.add_argument("-pn_dl", help="The number of times the median variability is required to trigger a detection.", type=float)
parser.add_argument("-m_dl", help="The number of times the median variability is required to trigger a detection.", type=float)
parser.add_argument("-bs", "--box-size", dest="bs", help="Size of the detection box in pixel^2.\nDefault: 5", default=5, nargs='?', type=int)
parser.add_argument("-tw", "--time-window", dest="tw", help="The duration of the time windows.\n Default: 100", default=100.0, nargs='?', type=float)
parser.add_argument("-inst", "--instrument", dest="inst", help="Type of detector", default='PN', nargs='?', type=str)
parser.add_argument("-gtr", "--good-time-ratio", dest="gtr", help="Ratio of acceptability for a time window. Shall be between 0.0 and 1.0.\nDefault: 1.0", default=1.0, nargs='?', type=float)

args = parser.parse_args()


pn_path = args.path + '/{}_{}_3_{}_PN/ds9_variable_sources.reg'.format(int(args.pn_dl), int(args.tw), args.gtr)
m_path = args.path + '/{}_{}_{}_{}_{}/ds9_variable_sources.reg'.format(int(args.m_dl), int(args.tw), args.bs, args.gtr, args.inst)

src_PN = Table(names=('ID', 'INST', 'RA', 'DEC') , dtype=('i2', 'U25', 'f8', 'f8'))
src_M = Table(names=('ID', 'INST', 'RA', 'DEC') , dtype=('i2', 'U25', 'f8', 'f8'))

# Implementing correlation table
corr_table = Table(names=('ID_1', 'INST_1', 'RA_1', 'DEC_1','ID_2', 'INST_2', 'RA_2', 'DEC_2', 'SEP') , dtype=('i2', 'U25', 'f8', 'f8','i2', 'U25', 'f8', 'f8', 'f8'))


pn_matches = 0
with open(pn_path,'r') as file:
    for line in file:
        if re.search("circle", line):
            pn_matches = pn_matches + 1
            #print(line)
            split = re.split('e |, ', line)
            src_PN.add_row((pn_matches, "PN", split[1], split[2]))


print("PN detections= ",pn_matches)

m_matches = 0
with open(m_path, 'r') as file:
    for line in file:
        if re.search("circle", line):
            m_matches = m_matches + 1
            #print(line)
            split = re.split('e |, ', line)
            src_M.add_row((m_matches, args.inst, split[1], split[2]))


print("M detections = ",m_matches)

# Checking correlation for the 3 EPIC
seperation_cut_off = 14
corr_table = check_correlation(src_PN, src_M, corr_table, sep_cutoff=seperation_cut_off)

# Sorting the table
corr_PN_M1 = corr_table[np.where((corr_table['INST_1'] == 'PN') & (corr_table['INST_2'] == args.inst))]

# Printing results
if len(corr_PN_M1) != 0:
    print("\n Correlation between EPIC-pn and EPIC-MOS1 \n")
    print(corr_PN_M1)

file_name = "match_count" + args.inst + ".txt"
match_count_file = open(file_name, 'w')
match_count_file.write(str(len(corr_PN_M1)))
match_count_file.close()

