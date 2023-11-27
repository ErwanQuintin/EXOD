"""
Eliminate any sources which lie within 1 arcmin of the bright source.
The list of bright sources from the XMM catalogue can be found in
bright_sources.csv table.
"""

from math import *
import argparse
import csv
import re

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord

from logger import logger


########################################################################
#                                                                      #
#          Parsing arguments                                           #
#                                                                      #
########################################################################
#
parser = argparse.ArgumentParser()
#
parser.add_argument(
    "-path", help="Path to the folder containing the XMM catalog", type=str
)

args = parser.parse_args()

fits_image_filename = "bright_sources.csv"
fits_image_filename = args.path + "/" + fits_image_filename


table_obs_id_list = pd.read_csv(fits_image_filename, usecols=["OBS_ID"])
table_ra_list = pd.read_csv(fits_image_filename, usecols=["SC_RA"])
table_dec_list = pd.read_csv(fits_image_filename, usecols=["SC_DEC"])
table_obs_id = table_obs_id_list["OBS_ID"]
table_ra = table_ra_list["SC_RA"]
table_dec = table_dec_list["SC_DEC"]

# Coulmn numbers in the triple_match.csv file

col_no_obs_id = 0
col_no_ra = 14
col_no_dec = 15


bright_src_distance_cutoff = 60  # 60 arcsec


with open("triple_match_obs.csv", "r") as csvinput:
    with open("triple_match_non_bright.csv", "w") as csvoutput:
        writer = csv.writer(csvoutput, lineterminator="\n")
        reader = csv.reader(csvinput)

        all = []
        row = next(reader)
        all.append(row)

        row_no = 0
        for row in reader:
            exod_ra = row[col_no_ra]
            exod_dec = row[col_no_dec]
            obs_id_ = row[col_no_obs_id]

            close_to_bright_src = 0
            iterator = 0
            for it in table_ra:
                if obs_id_ != table_obs_id[iterator]:
                    iterator = iterator + 1
                    continue
                else:
                    c1 = SkyCoord(
                        float(exod_ra), float(exod_dec), frame="fk5", unit="deg"
                    )
                    c2 = SkyCoord(
                        float(table_ra[iterator]),
                        float(table_dec[iterator]),
                        frame="fk5",
                        unit="deg",
                    )
                    sep = c1.separation(c2)
                    if sep.arcsecond < bright_src_distance_cutoff:
                        close_to_bright_src = 1
                        break
                    iterator = iterator + 1

            if close_to_bright_src == 0:
                all.append(row)
            row_no = row_no + 1

        writer.writerows(all)
