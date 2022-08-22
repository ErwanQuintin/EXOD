#!/usr/bin/env python3
# coding=utf-8


########################################################################
#                                                                      #
# EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System           #
#                                                                      #
# Generate bright_sources.csv table                                    #
#                                                                      #
# Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu                 #
#                                                                      #
########################################################################
"""
The following code generates a list of bright sources from the
XMM catalogue. The generated bright_sources.csv table
is used by remove_bright_sources.py. 
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


output_filename = "bright_sources.csv"
output_filename= args.path + "/" + output_filename


hdul = fits.open(fits_image_filename)
data = hdul[1].data
table_det_id = data.field('DETID')
table_src_id = data.field('SRCID')
table_obs_id = data.field('OBS_ID')
table_flux = data.field('EP_8_FLUX')
table_flux_err = data.field('EP_8_FLUX_ERR')
table_ra = data.field('SC_RA')
table_dec = data.field('SC_DEC')
table_pos_err = data.field('SC_POSERR')




with open(output_filename, 'w') as csvoutput:
    writer = csv.writer(csvoutput, lineterminator='\n')

    all = []
    row = []
    row.append('DETID')
    row.append('SRCID')
    row.append('OBS_ID')
    row.append('EP_8_FLUX')
    row.append('EP_8_FLUX_ERR')
    row.append('SC_RA')
    row.append('SC_DEC')
    row.append('SC_POSERR')
    all.append(row)

    it=0
    for detid in table_det_id:
        det_id = str(table_det_id[it])
        src_id = str(table_src_id[it])
        obs_id = str(table_obs_id[it])
        flux = table_flux[it]
        flux_err = str(table_flux_err[it])
        ra =str(table_ra[it])
        dec =str(table_dec[it])
        pos_err = str(table_pos_err[it])
        if(flux > 1e-11):
            row = []
            row.append(det_id)
            row.append(src_id)
            row.append("XMM " + obs_id)
            row.append(flux)
            row.append(flux_err)
            row.append(ra)
            row.append(dec)
            row.append(pos_err)
            all.append(row)
        it = it+1
    writer.writerows(all)

