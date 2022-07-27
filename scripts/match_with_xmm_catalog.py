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
table_obs_id = data.field('OBS_ID')
table_src_id = data.field('SRCID')
table_ra = data.field('SC_RA')
table_dec = data.field('SC_DEC')
table_flux = data.field('EP_8_FLUX')
table_var_flag = data.field('VAR_FLAG')
src_distance_cutoff = 30 # 30 arcsec


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
        row.append('XMM_RA')
        row.append('XMM_DEC')
        row.append('EXOD_XMM_SEP')
        row.append('XMM_SRCID')
        row.append('XMM_EP8_FLUX')
        row.append('XMM_PIPELINE_VARIABLE')
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
            if(row[col_no_type] == "T"):
                center_ra = str( ((float(pn_ra) + float(m1_ra) + float(m2_ra) )/3.0))
                center_dec = str(((float(pn_dec) + float(m1_dec) + float(m2_dec)) / 3.0))
            if(row[col_no_type] == "D"):
                if(m2_ra =="-"):
                    center_ra = str( ((float(pn_ra) + float(m1_ra)  )/2.0))
                    center_dec = str(((float(pn_dec) + float(m1_dec) ) / 2.0))
                if(m1_ra =="-"):
                    center_ra = str( ((float(pn_ra) + float(m2_ra)  )/2.0))
                    center_dec = str(((float(pn_dec) + float(m2_dec) ) / 2.0))
                if(pn_ra =="-"):
                    center_ra = str( ((float(m1_ra) + float(m2_ra)  )/2.0))
                    center_dec = str(((float(m1_dec) + float(m2_dec) ) / 2.0))
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


            split_row = re.split(' ', row[0])
            row.append(center_ra)
            row.append(center_dec)
            iterator =0
            closest_sep=-1
            closest_ra="NA"
            closest_dec="NA"
            closest_src_id ="NA"
            closest_flux ="NA"
            closest_var_flag ="No_match"
            print("working on OBS ID", split_row[1])
            for it in table_obs_id:
                if (it != split_row[1]):
                    iterator = iterator + 1
                    continue
                else:
                    c1 = SkyCoord(float(center_ra), float(center_dec), frame='fk5', unit='deg')
                    c2 = SkyCoord(float(table_ra[iterator]), float(table_dec[iterator]), frame='fk5', unit='deg')
                    sep = c1.separation(c2)
                    if(sep.arcsecond <= src_distance_cutoff):
                        if(closest_sep==-1):
                            closest_sep=(round(sep.arcsecond,2))
                            closest_ra=float(table_ra[iterator])
                            closest_dec=float(table_dec[iterator])
                            closest_src_id = table_src_id[iterator]
                            closest_flux = table_flux[iterator]
                            if(table_var_flag[iterator] == True): 
                                closest_var_flag ="Yes_and_variable"
                            else:
                                closest_var_flag ="Yes_and_not_variable"
                        else:
                            if(closest_sep > sep.arcsecond):
                                closest_sep=(round(sep.arcsecond,2))
                                closest_ra=-float(table_ra[iterator])
                                closest_dec=float(table_dec[iterator])
                                closest_src_id = table_src_id[iterator]
                                closest_flux = table_flux[iterator]
                                if(table_var_flag[iterator] == True): 
                                    closest_var_flag ="Yes_and_variable"
                                else:
                                    closest_var_flag ="Yes_and_not_variable"
                    iterator = iterator + 1

            if(closest_sep!=-1):
                row.append(closest_ra)
                row.append(closest_dec)
                row.append(str(round(closest_sep,2)))
                row.append("XMM" +str(closest_src_id))
                row.append(closest_flux)
                row.append(closest_var_flag)
            else:
                row.append(closest_ra)
                row.append(closest_dec)
                row.append(str(round(closest_sep,2)))
                row.append("XMM" +str(closest_src_id))
                row.append(closest_flux)
                row.append(closest_var_flag)
            all.append(row)

        writer.writerows(all)

