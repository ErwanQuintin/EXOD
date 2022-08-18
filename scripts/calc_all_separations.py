#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Cross correlate variable sources                                      #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################
"""
Cross correlate variable sources across the 3 EPIC cameras and
determine distances between the sources and write results into
separation_file.csv which will be used by detect_final_corr.py

"""

import re
import sys
import os
import numpy as np
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table
import argparse


seperation_cut_off = 600
########################################################################
#                                                                      #
# Determine the seperation between two sources and label them as       #
# correlated if separation is less than seperation_cut_off            #
#                                                                      #
########################################################################

def check_correlation(src_1, src_2, corr_tab):
    """
	Function checking the correlation sources between two source lists

	@param  src_1:      The source list of first detector
	@param  src_2:      The source list of second detector
	@param  corr_table: The correlation table

	@return: The correlation table appended
	"""

    for i in range(len(src_1)):
        for j in range(len(src_2)):
            c1 = SkyCoord(src_1['RA'][i], src_1['DEC'][i], frame='fk5', unit='deg')
            c2 = SkyCoord(src_2['RA'][j], src_2['DEC'][j], frame='fk5', unit='deg')
            sep = c1.separation(c2)
            if sep.arcsecond < seperation_cut_off:
                corr_tab.add_row([src_1['ID'][i], src_1['INST'][i], src_1['RA'][i], src_1['DEC'][i], src_2['ID'][j], src_2['INST'][j], src_2['RA'][j], src_2['DEC'][j], round(sep.arcsecond,2)])

    return corr_tab


########################################################################
#                                                                      #
#          Parsing arguments                                           #
#                                                                      #
########################################################################
##
parser = argparse.ArgumentParser()

# Variability parameters
parser.add_argument("-path", help="Path to the folder containing the observation files", type=str)
parser.add_argument("-pn_dl", help="The number of times the median variability is required to trigger a detection.", type=float)
parser.add_argument("-m1_dl", help="The number of times the median variability is required to trigger a detection.", type=float)
parser.add_argument("-m2_dl", help="The number of times the median variability is required to trigger a detection.", type=float)
parser.add_argument("-bs", "--box-size", dest="bs", help="Size of the detection box in pixel^2.\nDefault: 5", default=5, nargs='?', type=int)
parser.add_argument("-tw", "--time-window", dest="tw", help="The duration of the time windows.\n Default: 100", default=100.0, nargs='?', type=float)
parser.add_argument("-inst", "--instrument", dest="inst", help="Type of detector", default='PN', nargs='?', type=str)
parser.add_argument("-obs", "--observation", dest="obs", help="Observation ID", default=None, nargs='?', type=str)
parser.add_argument("-gtr", "--good-time-ratio", dest="gtr", help="Ratio of acceptability for a time window. Shall be between 0.0 and 1.0.\nDefault: 1.0", default=1.0, nargs='?', type=float)

args = parser.parse_args()

print("pn_dl =", args.pn_dl," m1_dl =", args.m1_dl," m2_dl =", args.m2_dl)


pn_path = args.path + '/{}_{}_3_{}_PN/ds9_variable_sources.reg'.format(int(args.pn_dl), int(args.tw), args.gtr)
mn1_path = args.path + '/{}_{}_{}_{}_M1/ds9_variable_sources.reg'.format(int(args.m1_dl), int(args.tw), args.bs, args.gtr)
mn2_path = args.path + '/{}_{}_{}_{}_M2/ds9_variable_sources.reg'.format(int(args.m2_dl), int(args.tw), args.bs, args.gtr)


src_PN = Table(names=('ID', 'INST', 'RA', 'DEC') , dtype=('i2', 'U25', 'f8', 'f8'))
src_M1 = Table(names=('ID', 'INST', 'RA', 'DEC') , dtype=('i2', 'U25', 'f8', 'f8'))
src_M2 = Table(names=('ID', 'INST', 'RA', 'DEC') , dtype=('i2', 'U25', 'f8', 'f8'))

# Implementing correlation table
corr_table = Table(names=('ID_1', 'INST_1', 'RA_1', 'DEC_1','ID_2', 'INST_2', 'RA_2', 'DEC_2', 'SEP') , dtype=('i2', 'U25', 'f8', 'f8','i2', 'U25', 'f8', 'f8', 'f8'))


pn_matches = 0
try:
    with open(pn_path, 'r') as file:
        for line in file:
            if re.search("circle", line):
                pn_matches = pn_matches + 1
                # print(line)
                split = re.split('e |, ', line)
                src_PN.add_row((pn_matches, "PN", split[1], split[2]))
except:
    print("WARNING:: Dint find file:", pn_path)
print("PN detections = ",pn_matches)

m1_matches = 0
try:
    with open(mn1_path, 'r') as file:
        for line in file:
            if re.search("circle", line):
                m1_matches = m1_matches + 1
                # print(line)
                split = re.split('e |, ', line)
                src_M1.add_row((m1_matches, "M1", split[1], split[2]))
except:
    print("WARNING:: Dint find file:", mn1_path)

print("M1 detections = ",m1_matches)

m2_matches = 0
try:
    with open(mn2_path, 'r') as file:
        for line in file:
            if re.search("circle", line):
                m2_matches = m2_matches + 1
                # print(line)
                split = re.split('e |, ', line)
                src_M2.add_row((m2_matches, "M2", split[1], split[2]))
except:
    print("WARNING:: Dint find file:", mn2_path)

print("M2 detections = ",m2_matches)
#print(src_M2)




# Checking correlation for the 3 EPIC
corr_table = check_correlation(src_PN, src_M1, corr_table)
corr_table = check_correlation(src_M1, src_M2, corr_table)
corr_table = check_correlation(src_PN, src_M2, corr_table)

# Sorting the table
corr_PN_M1 = corr_table[np.where((corr_table['INST_1'] == 'PN') & (corr_table['INST_2'] == 'M1'))]
corr_M1_M2 = corr_table[np.where((corr_table['INST_1'] == 'M1') & (corr_table['INST_2'] == 'M2'))]
corr_PN_M2 = corr_table[np.where((corr_table['INST_1'] == 'PN') & (corr_table['INST_2'] == 'M2'))]

# Printing results
if len(corr_PN_M1) != 0:
    print("\n Correlation between EPIC-pn and EPIC-MOS1 \n")
    corr_PN_M1.sort('SEP')
    print(corr_PN_M1)

if len(corr_M1_M2) != 0:
    print("\n Correlation between EPIC-MOS1 and EPIC-MOS2 \n")
    corr_M1_M2.sort('SEP')
    print(corr_M1_M2)

if len(corr_PN_M2) != 0:
    print("\n Correlation between EPIC-pn and EPIC-MOS2 \n")
    corr_PN_M2.sort('SEP')
    print(corr_PN_M2)



if (os.path.exists('./separation_file.csv')):
        master_file = open("separation_file.csv", 'a')
else:
        master_file = open("separation_file.csv", 'w')
        master_file.write("OBS_ID,class,PN_REGION_NUMBER,PN_RA,PN_DEC,M1_REGION_NUMBER,M1_RA,M1_DEC,M2_REGION_NUMBER,M2_RA, M2_DEC,PN_M1_SEP,PN_M2_SEP, M1_M2_SEP\n")


if(True):
    for ii in range(len(corr_PN_M1)):
        line = """XMM {0},D,{1},{2},{3},{4},{5},{6},-,-,-,{7},-,-\n""".format(str(args.obs),corr_PN_M1[ii]['ID_1'],
                                                                                     corr_PN_M1[ii]['RA_1'],
                                                                                     corr_PN_M1[ii]['DEC_1'],
                                                                                     corr_PN_M1[ii]['ID_2'],
                                                                                     corr_PN_M1[ii]['RA_2'],
                                                                                     corr_PN_M1[ii]['DEC_2'],
                                                                                     corr_PN_M1[ii]['SEP'])
        master_file.write(line)

    for ii in range(len(corr_PN_M2)):
        line = """XMM {0},D,{1},{2},{3},-,-,-,{4},{5},{6},-,{7},-\n""".format(str(args.obs),corr_PN_M2[ii]['ID_1'],
                                                                                     corr_PN_M2[ii]['RA_1'],
                                                                                     corr_PN_M2[ii]['DEC_1'],
                                                                                     corr_PN_M2[ii]['ID_2'],
                                                                                     corr_PN_M2[ii]['RA_2'],
                                                                                     corr_PN_M2[ii]['DEC_2'],
                                                                                     corr_PN_M2[ii]['SEP'])
        master_file.write(line)

    for ii in range(len(corr_M1_M2)):
        line = """XMM {0},D,-,-,-,{1},{2},{3},{4},{5},{6},-,-,{7}\n""".format(str(args.obs),corr_M1_M2[ii]['ID_1'],
                                                                                     corr_M1_M2[ii]['RA_1'],
                                                                                     corr_M1_M2[ii]['DEC_1'],
                                                                                     corr_M1_M2[ii]['ID_2'],
                                                                                     corr_M1_M2[ii]['RA_2'],
                                                                                     corr_M1_M2[ii]['DEC_2'],
                                                                                     corr_M1_M2[ii]['SEP'])
        master_file.write(line)

master_file.close()

