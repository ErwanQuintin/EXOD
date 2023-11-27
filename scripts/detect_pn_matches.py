#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Independetly detect the number of PN matches                         #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################
"""
Independetly detect the number of PN matches from a given past EXOD run

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

# Third-party imports
import argparse


########################################################################
#                                                                      #
#          Parsing arguments                                           #
#                                                                      #
########################################################################
#
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
    "-gtr",
    "--good-time-ratio",
    dest="gtr",
    help="Ratio of acceptability for a time window. Shall be between 0.0 and 1.0.\nDefault: 1.0",
    default=1.0,
    nargs="?",
    type=float,
)

args = parser.parse_args()


pn_path = args.path + "/{}_{}_3_{}_PN/ds9_variable_sources.reg".format(
    int(args.pn_dl), int(args.tw), args.gtr
)

# pn_path = args.path + "\/" + '{}_{}_3_{}_PN\ds9_variable_sources.reg'.format(int(args.pn_dl), int(args.tw), args.gtr)

src_PN = Table(names=("ID", "INST", "RA", "DEC"), dtype=("i2", "U25", "f8", "f8"))


pn_matches = 0
with open(pn_path, "r") as file:
    for line in file:
        if re.search("circle", line):
            pn_matches = pn_matches + 1
            # print(line)
            split = re.split("e |, ", line)
            src_PN.add_row((pn_matches, "PN", split[1], split[2]))


print("PN matches = ", pn_matches)

match_count_file = open("match_count_pn.txt", "w")
match_count_file.write(str(pn_matches))
match_count_file.close()
