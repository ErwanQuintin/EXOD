#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Cross match sources with 4XMM_DR11cat_v1 catalog                     #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################
"""
 Cross match sources with 4XMM_DR11cat_v1 catalog
 If it's a triple or double, the centroid is matched with the closest XMM source
 If it's a single source, it gets matched conventionally.
"""
from math import *
import csv
import pandas as pd
import math
import numpy as np
from astropy.io import fits
import re
from astropy.coordinates import SkyCoord
import argparse

########################################################################
#                                                                      #
#          Parsing arguments                                           #
#                                                                      #
########################################################################
#
parser = argparse.ArgumentParser()
#
parser.add_argument("-path", help="Path to the folder containing the XMM catalog", type=str)

args = parser.parse_args()

fits_image_filename = "bright_sources.fits"
fits_image_filename= args.path + "/" + fits_image_filename
hdul = fits.open(fits_image_filename)
data = hdul[1].data
table_obs_id = data.field('OBS_ID')
table_ra = data.field('SC_RA')
table_dec = data.field('SC_DEC')

#Coulmn numbers in the triple_match.csv file

col_no_obs_id=0
col_no_ra=14
col_no_dec=15


bright_src_distance_cutoff = 60 # 60 arcsec


with open('triple_match_obs.csv','r') as csvinput:
    with open('triple_match_non_bright.csv', 'w') as csvoutput:
        writer = csv.writer(csvoutput, lineterminator='\n')
        reader = csv.reader(csvinput)

        all = []
        row = next(reader)
        all.append(row)

        row_no=0
        for row in reader:
            exod_ra= row[col_no_ra]
            exod_dec= row[col_no_dec]
            split_row = re.split(' ', row[col_no_obs_id])
            obs_id_ =  split_row[1]

            close_to_bright_src=0
            iterator=0
            for it in table_ra:
                if(obs_id_ != table_obs_id[iterator]):
                    iterator = iterator + 1
                    continue
                else:
                    c1 = SkyCoord(float(exod_ra), float(exod_dec), frame='fk5', unit='deg')
                    c2 = SkyCoord(float(table_ra[iterator]), float(table_dec[iterator]), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    if(sep.arcsecond < bright_src_distance_cutoff):
                        close_to_bright_src=1
                        break
                    iterator = iterator + 1

            if(close_to_bright_src==0):
                all.append(row)
            row_no=row_no+1

        writer.writerows(all)

