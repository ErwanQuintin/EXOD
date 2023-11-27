#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Collect data for independent light curve generation                  #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################
"""
Collect data for independent light curve generation. 
This allows us to decouple the light curve generation 
from the variable source dection.
"""

import re
import sys
import os
import numpy as np
import argparse


########################################################################
#                                                                      #
#          Parsing arguments                                           #
#                                                                      #
########################################################################
##
parser = argparse.ArgumentParser()

# Variability parameters
parser.add_argument(
    "-path", help="Path to the folder containing the observation files", type=str
)
parser.add_argument(
    "-pn_dl",
    help="The number of times the median variability is required to trigger a detection.",
    type=float,
)
parser.add_argument(
    "-m1_dl",
    help="The number of times the median variability is required to trigger a detection.",
    type=float,
)
parser.add_argument(
    "-m2_dl",
    help="The number of times the median variability is required to trigger a detection.",
    type=float,
)
parser.add_argument(
    "-bs",
    "--box-size",
    dest="bs",
    help="Size of the detection box in pixel^2.\nDefault: 5",
    default=5,
    nargs="?",
    type=int,
)
parser.add_argument(
    "-tw",
    "--time-window",
    dest="tw",
    help="The duration of the time windows.\n Default: 100",
    default=100.0,
    nargs="?",
    type=float,
)
parser.add_argument(
    "-inst",
    "--instrument",
    dest="inst",
    help="Type of detector",
    default="PN",
    nargs="?",
    type=str,
)
parser.add_argument(
    "-obs",
    "--observation",
    dest="obs",
    help="Observation ID",
    default=None,
    nargs="?",
    type=str,
)
parser.add_argument(
    "-gtr",
    "--good-time-ratio",
    dest="gtr",
    help="Ratio of acceptability for a time window. Shall be between 0.0 and 1.0.\nDefault: 1.0",
    default=1.0,
    nargs="?",
    type=float,
)

args = parser.parse_args()

print("pn_dl =", args.pn_dl, " m1_dl =", args.m1_dl, " m2_dl =", args.m2_dl)

pn_path = args.path + "/{}_{}_3_{}_PN/ds9_variable_sources.reg".format(
    int(args.pn_dl), int(args.tw), args.gtr
)
mn1_path = args.path + "/{}_{}_{}_{}_M1/ds9_variable_sources.reg".format(
    int(args.m1_dl), int(args.tw), args.bs, args.gtr
)
mn2_path = args.path + "/{}_{}_{}_{}_M2/ds9_variable_sources.reg".format(
    int(args.m2_dl), int(args.tw), args.bs, args.gtr
)

pn_matches = 0
try:
    with open(pn_path, "r") as file:
        for line in file:
            if re.search("circle", line):
                pn_matches = pn_matches + 1
except:
    print("WARNING:: Dint find file:", pn_path)
print("PN detections = ", pn_matches)

m1_matches = 0
try:
    with open(mn1_path, "r") as file:
        for line in file:
            if re.search("circle", line):
                m1_matches = m1_matches + 1
except:
    print("WARNING:: Dint find file:", mn1_path)

print("M1 detections = ", m1_matches)

m2_matches = 0
try:
    with open(mn2_path, "r") as file:
        for line in file:
            if re.search("circle", line):
                m2_matches = m2_matches + 1
except:
    print("WARNING:: Dint find file:", mn2_path)

print("M2 detections = ", m2_matches)


file_name = "match_countPN.txt"
match_count_file = open(file_name, "w")
match_count_file.write(str(pn_matches + 1))
match_count_file.close()

file_name = "match_countM1.txt"
match_count_file = open(file_name, "w")
match_count_file.write(str(m1_matches + 1))
match_count_file.close()


file_name = "match_countM2.txt"
match_count_file = open(file_name, "w")
match_count_file.write(str(m2_matches + 1))
match_count_file.close()
