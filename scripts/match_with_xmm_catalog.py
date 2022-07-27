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

fits_image_filename = "4XMM_DR11cat_v1.0.fits"
fits_image_filename= args.path + "/" + fits_image_filename
hdul = fits.open(fits_image_filename)
data = hdul[1].data
obs_id = data.field('OBS_ID')
table_ra = data.field('SC_RA')
table_dec = data.field('SC_DEC')

#Coulmn numbers in the triple_match.csv file

col_no_obs_id=0
col_no_type=1
col_no_pn_ra=3
col_no_pn_dec=4
col_no_m1_ra=6
col_no_m1_dec=7
col_no_m2_ra=9
col_no_m2_dec=10
col_no_pn_m1=11
col_no_pn_m2=12
col_no_m1_m2=13


with open('triple_match.csv','r') as csvinput:
    with open('triple_match_obs.csv', 'w') as csvoutput:
        writer = csv.writer(csvoutput, lineterminator='\n')
        reader = csv.reader(csvinput)

        all = []
        row = next(reader)
        row.append('EXOD_RA')
        row.append('EXOD_DEC')
        row.append('EXOD_RA_MEDIAN')
        row.append('EXOD_DEC_MEDIAN')
        row.append('EXOD_PN_SEP')
        row.append('EXOD_M1_SEP')
        row.append('EXOD_M2_SEP')
        row.append('EXOD_PN_SEP_MED')
        row.append('EXOD_M1_SEP_MED')
        row.append('EXOD_M2_SEP_MED')
        row.append('XMM_RA')
        row.append('XMM_DEC')
        row.append('EXOD_XMM_SEP')
        all.append(row)

        for row in reader:
            pn_ra= row[col_no_pn_ra]
            m1_ra= row[col_no_m1_ra]
            m2_ra= row[col_no_m2_ra]
            pn_dec= row[col_no_pn_dec]
            m1_dec= row[col_no_m1_dec]
            m2_dec= row[col_no_m2_dec]
            obj_class= row[col_no_type]
            center_ra = '-'
            center_dec = '-'
            center_pn_sep = '-'
            center_m1_sep = '-'
            center_m2_sep = '-'
            center_ra_med = '-'
            center_dec_med = '-'
            center_pn_sep_med = '-'
            center_m1_sep_med = '-'
            center_m2_sep_med = '-'
            if(row[col_no_type] == "T"):
                center_ra = str( ((float(pn_ra) + float(m1_ra) + float(m2_ra) )/3.0))
                center_dec = str(((float(pn_dec) + float(m1_dec) + float(m2_dec)) / 3.0))
                c1 = SkyCoord(float(center_ra), float(center_dec), frame='fk5', unit='deg')
                c2 = SkyCoord(float(pn_ra), float(pn_dec), frame='fk5', unit='deg')
                sep = c1.separation(c2)
                center_pn_sep = str(round(sep.arcsecond, 2))
                c2 = SkyCoord(float(m1_ra), float(m1_dec), frame='fk5', unit='deg')
                sep = c1.separation(c2)
                center_m1_sep = str(round(sep.arcsecond, 2))
                c2 = SkyCoord(float(m2_ra), float(m2_dec), frame='fk5', unit='deg')
                sep = c1.separation(c2)
                center_m2_sep = str(round(sep.arcsecond, 2))

                center_ra_med = str( (np.median([float(pn_ra) , float(m1_ra) , float(m2_ra)]) ))
                center_dec_med =  str( (np.median([float(pn_dec) , float(m1_dec) , float(m2_dec)]) ))
                c1 = SkyCoord(float(center_ra_med), float(center_dec_med), frame='fk5', unit='deg')
                c2 = SkyCoord(float(pn_ra), float(pn_dec), frame='fk5', unit='deg')
                sep = c1.separation(c2)
                center_pn_sep_med = str(round(sep.arcsecond, 2))
                c2 = SkyCoord(float(m1_ra), float(m1_dec), frame='fk5', unit='deg')
                sep = c1.separation(c2)
                center_m1_sep_med = str(round(sep.arcsecond, 2))
                c2 = SkyCoord(float(m2_ra), float(m2_dec), frame='fk5', unit='deg')
                sep = c1.separation(c2)
                center_m2_sep_med = str(round(sep.arcsecond, 2))

            if(row[col_no_type] == "D"):
                if(m2_ra =="-"):
                    center_ra = str( ((float(pn_ra) + float(m1_ra)  )/2.0))
                    center_dec = str(((float(pn_dec) + float(m1_dec) ) / 2.0))
                    c1 = SkyCoord(float(center_ra), float(center_dec), frame='fk5', unit='deg')
                    c2 = SkyCoord(float(pn_ra), float(pn_dec), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    center_pn_sep = str(round(sep.arcsecond, 2))
                    c2 = SkyCoord(float(m1_ra), float(m1_dec), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    center_m1_sep = str(round(sep.arcsecond, 2))
                if(m1_ra =="-"):
                    center_ra = str( ((float(pn_ra) + float(m2_ra)  )/2.0))
                    center_dec = str(((float(pn_dec) + float(m2_dec) ) / 2.0))
                    c1 = SkyCoord(float(center_ra), float(center_dec), frame='fk5', unit='deg')
                    c2 = SkyCoord(float(pn_ra), float(pn_dec), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    center_pn_sep = str(round(sep.arcsecond, 2))
                    c2 = SkyCoord(float(m2_ra), float(m2_dec), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    center_m2_sep = str(round(sep.arcsecond, 2))
                if(pn_ra =="-"):
                    center_ra = str( ((float(m1_ra) + float(m2_ra)  )/2.0))
                    center_dec = str(((float(m1_dec) + float(m2_dec) ) / 2.0))
                    c1 = SkyCoord(float(center_ra), float(center_dec), frame='fk5', unit='deg')
                    c2 = SkyCoord(float(m1_ra), float(m1_dec), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    center_m1_sep = str(round(sep.arcsecond, 2))
                    c2 = SkyCoord(float(m2_ra), float(m2_dec), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    center_m2_sep = str(round(sep.arcsecond, 2))
                center_ra_med = center_ra
                center_dec_med = center_dec
                center_pn_sep_med = center_pn_sep
                center_m1_sep_med = center_m1_sep
                center_m2_sep_med = center_m2_sep
            if(row[col_no_type] == "S"):
                if(m2_ra !="-"):
                    center_ra = m2_ra
                    center_dec = m2_dec
                if(m1_ra !="-"):
                    center_ra = m1_ra
                    center_dec = m1_dec
                if(pn_ra !="-"):
                    center_ra = pn_ra
                    center_dec = pn_dec
                center_ra_med = center_ra
                center_dec_med = center_dec
                center_pn_sep_med = center_pn_sep
                center_m1_sep_med = center_m1_sep
                center_m2_sep_med = center_m2_sep


            split_row = re.split(' ', row[0])
            row.append(center_ra)
            row.append(center_dec)
            row.append(center_ra_med)
            row.append(center_dec_med)
            row.append(center_pn_sep)
            row.append(center_m1_sep)
            row.append(center_m2_sep)
            row.append(center_pn_sep_med)
            row.append(center_m1_sep_med)
            row.append(center_m2_sep_med)
            jj =0
            for it in obs_id:
                if(it == split_row[1]):
                    row.append(table_ra[jj])
                    row.append(table_dec[jj])
                    c1 = SkyCoord(float(center_ra), float(center_dec), frame='fk5', unit='deg')
                    c2 = SkyCoord(float(table_ra[jj]), float(table_dec[jj]), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    row.append( str(round(sep.arcsecond,2)))
                    all.append(row)
                    break
                jj=jj+1

        writer.writerows(all)

