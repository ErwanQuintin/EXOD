#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
#  Determine triple/double correlated sources                           #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################
"""
Use the cross correlate variable sources across the 3 EPIC cameras
from the previous step and organise the results into triple matches
and double matches. Sources which are not cross correlated are
preseved and labled as singles.
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

from detect_corr import check_correlation


seperation_cut_off = 25
include_doubles=True
include_singles=True

########################################################################
#                                                                      #
#  Take 3 sources across the 3 EPIC cameras and check if               #
#   they are they same source. If so lael it as triple match           #
#                                                                      #
########################################################################


def check_triple(corr_1, corr_2, corr_3):
    """
	Function checking correlations between 3 correlation tables

	@param  corr_1: The first correlation table
	@param  corr_2: The second correlation table
	@param  corr_3: The third correlation table

	@return: A list with ID and INST of triple correlation sources
	"""

    src_cand = []
    triple = []

    for i in range(len(corr_1)):
        src_cand.append([corr_1['ID_1'][i], corr_1['INST_1'][i], corr_1['ID_2'][i], corr_1['INST_2'][i]])


    for cand in src_cand:
        t1 = corr_3[np.where((corr_3['ID_1'] == cand[0]) & (corr_3['INST_1'] == cand[1]))]
        t2 = corr_2[np.where((corr_2['ID_1'] == cand[2]) & (corr_2['INST_1'] == cand[3]))]

        for line1 in t1:
            for line2 in t2:
                if line1['ID_2'] == line2['ID_2'] and line1['INST_2'] == line2['INST_2']:
                    triple.append(
                        [line1['ID_1'], line1['INST_1'], line1['RA_1'], line1['DEC_1'], line2['ID_1'], line2['INST_1'], line2['RA_1'], line2['DEC_1'], line1['ID_2'], line1['INST_2'], line1['RA_2'], line1['DEC_2']])

    return triple


########################################################################
#                                                                      #
#          Parsing arguments                                           #
#                                                                      #
########################################################################
#
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

# Checking for triple correlation
src_triple = []
if len(corr_PN_M1) != 0 and len(corr_M1_M2) != 0 and len(corr_PN_M2) != 0:
    src_triple = check_triple(corr_PN_M1, corr_M1_M2, corr_PN_M2)


if len(src_triple) !=0:
    print("\n Correlation between the 3 EPIC detectors \n")
    print(src_triple)

if (os.path.exists('./triple_match.csv')):
        master_file = open("triple_match.csv", 'a')
else:
        master_file = open("triple_match.csv", 'w')
        master_file.write("OBS_ID,class,PN_REGION_NUMBER,PN_RA,PN_DEC,M1_REGION_NUMBER,M1_RA,M1_DEC,M2_REGION_NUMBER,M2_RA, M2_DEC,PN_M1_SEP,PN_M2_SEP, M1_M2_SEP\n")

for it in range(len(src_triple)):
    pn_c = SkyCoord(src_triple[it][2],src_triple[it][3], frame='fk5', unit='deg')
    m1_c = SkyCoord(src_triple[it][6],src_triple[it][7], frame='fk5', unit='deg')
    m2_c = SkyCoord(src_triple[it][10],src_triple[it][11], frame='fk5', unit='deg')
    pn_m1_sep_d = pn_c.separation(m1_c)
    pn_m2_sep_d = pn_c.separation(m2_c)
    m1_m2_sep_d = m1_c.separation(m2_c)

    pn_m1_sep= round(pn_m1_sep_d.arcsecond, 2)
    pn_m2_sep= round(pn_m2_sep_d.arcsecond, 2)
    m1_m2_sep= round(m1_m2_sep_d.arcsecond, 2)


    line = """XMM {0},T,{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}\n""".format(str(args.obs), src_triple[it][0],src_triple[it][2],src_triple[it][3],src_triple[it][4],src_triple[it][6],src_triple[it][7],src_triple[it][8],src_triple[it][10],src_triple[it][11],pn_m1_sep,pn_m2_sep,m1_m2_sep)
    master_file.write(line)

if(include_doubles):
    for ii in range(len(corr_PN_M1)):
        present = 0;
        id = corr_PN_M1[ii]['ID_1']
        for jj in range(len(src_triple)):
            #print("looking for ID", id, "table ID ",  src_triple[jj][0])
            if(id ==  src_triple[jj][0]):
                present = present + 1
        if(present == 0):
            line = """XMM {0},D,{1},{2},{3},{4},{5},{6},-,-,-,{7},-,-\n""".format(str(args.obs),corr_PN_M1[ii]['ID_1'],
                                                                                     corr_PN_M1[ii]['RA_1'],
                                                                                     corr_PN_M1[ii]['DEC_1'],
                                                                                     corr_PN_M1[ii]['ID_2'],
                                                                                     corr_PN_M1[ii]['RA_2'],
                                                                                     corr_PN_M1[ii]['DEC_2'],
                                                                                     corr_PN_M1[ii]['SEP'])
            master_file.write(line)

    for ii in range(len(corr_PN_M2)):
        present = 0;
        id = corr_PN_M2[ii]['ID_1']
        for jj in range(len(src_triple)):
            #print("looking for ID", id, "table ID ",  src_triple[jj][0])
            if(id ==  src_triple[jj][0]):
                present = present + 1
        if(present == 0):
            line = """XMM {0},D,{1},{2},{3},-,-,-,{4},{5},{6},-,{7},-\n""".format(str(args.obs),corr_PN_M2[ii]['ID_1'],
                                                                                     corr_PN_M2[ii]['RA_1'],
                                                                                     corr_PN_M2[ii]['DEC_1'],
                                                                                     corr_PN_M2[ii]['ID_2'],
                                                                                     corr_PN_M2[ii]['RA_2'],
                                                                                     corr_PN_M2[ii]['DEC_2'],
                                                                                     corr_PN_M2[ii]['SEP'])
            master_file.write(line)

    for ii in range(len(corr_M1_M2)):
        present = 0;
        id = corr_M1_M2[ii]['ID_1']
        for jj in range(len(src_triple)):
            #print("looking for ID", id, "table ID ",  src_triple[jj][0])
            if(id ==  src_triple[jj][0]):
                present = present + 1
        if(present == 0):
            line = """XMM {0},D,-,-,-,{1},{2},{3},{4},{5},{6},-,-,{7}\n""".format(str(args.obs),corr_M1_M2[ii]['ID_1'],
                                                                                     corr_M1_M2[ii]['RA_1'],
                                                                                     corr_M1_M2[ii]['DEC_1'],
                                                                                     corr_M1_M2[ii]['ID_2'],
                                                                                     corr_M1_M2[ii]['RA_2'],
                                                                                     corr_M1_M2[ii]['DEC_2'],
                                                                                     corr_M1_M2[ii]['SEP'])
            master_file.write(line)


if(include_singles):

    #print(src_PN)
    for ii in range(len(src_PN)):
        present = 0;
        id = src_PN[ii]['ID']
        for jj in range(len(corr_PN_M1)):
            #print("looking for ID", id, "table ID ", corr_PN_M1[jj]['ID_1'])
            if (id == corr_PN_M1[jj]['ID_1']):
                present = present + 1
        for jj in range(len(corr_PN_M2)):
            #print("looking for ID", id, "table ID ", corr_PN_M2[jj]['ID_1'])
            if (id == corr_PN_M2[jj]['ID_1']):
                present = present + 1
        if (present == 0):
            line = """XMM {0},S,{1},{2},{3},-,-,-,-,-,-,-,-,-\n""".format(str(args.obs), src_PN[ii]['ID'],
                                                                    src_PN[ii]['RA'],
                                                                    src_PN[ii]['DEC'])
            master_file.write(line)

    #print(src_M1)
    for ii in range(len(src_M1)):
        present = 0;
        id = src_M1[ii]['ID']
        for jj in range(len(corr_M1_M2)):
            #print("looking for ID", id, "table ID ", corr_M1_M2[jj]['ID_1'])
            if (id == corr_M1_M2[jj]['ID_1']):
                present = present + 1
        for jj in range(len(corr_PN_M2)):
            #print("looking for ID", id, "table ID ", corr_PN_M1[jj]['ID_2'])
            if (id == corr_PN_M2[jj]['ID_2']):
                present = present + 1
        if (present == 0):
            line = """XMM {0},S,-,-,-,{1},{2},{3},-,-,-,-,-,-\n""".format(str(args.obs), src_M1[ii]['ID'],
                                                                            src_M1[ii]['RA'],
                                                                            src_M1[ii]['DEC'])
            master_file.write(line)

    #print(src_M2)
    for ii in range(len(src_M2)):
        present = 0;
        id = src_M2[ii]['ID']
        for jj in range(len(corr_M1_M2)):
            #print("looking for ID", id, "table ID ", corr_M1_M2[jj]['ID_2'])
            if (id == corr_M1_M2[jj]['ID_2']):
                present = present + 1
        for jj in range(len(corr_PN_M2)):
            #print("looking for ID", id, "table ID ", corr_PN_M2[jj]['ID_2'])
            if (id == corr_PN_M2[jj]['ID_2']):
                present = present + 1
        if (present == 0):
            line = """XMM {0},S,-,-,-,-,-,-,{1},{2},{3},-,-,-\n""".format(str(args.obs), src_M2[ii]['ID'],
                                                                            src_M2[ii]['RA'],
                                                                            src_M2[ii]['DEC'])
            master_file.write(line)

master_file.close()
