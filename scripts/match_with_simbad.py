"""
Query SIMBAD and find the closest 5 matches 
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System

Given a RA, DEC query SIMBAD and find the closest 5 matches that fall
within 30 arcsec of the detection and collate the results into a table.
"""
import os
import csv

import astropy.coordinates as coord
import astropy.units as u
from astroquery.simbad import Simbad
from astropy.io import ascii
from astropy.coordinates import SkyCoord

from logger import logger


########################################################################
#                                                                      #
# Given a RA, DEC query SIMBAD and find the closest 5 matches          #
# that fall within 30 arcsec of the detection                          #
#                                                                      #
########################################################################


def find_simbad(my_ra, my_dec):
    """
    Function to query SIMBAD
    @param  my_ra:  The RA of the object
    @param  my_dec: The DEC of the object
    @return: A list of 5 closest matches
    """

    match1_sep = "NA"
    match2_sep = "NA"
    match3_sep = "NA"
    match4_sep = "NA"
    match5_sep = "NA"

    match1_obj_name = "NA"
    match2_obj_name = "NA"
    match3_obj_name = "NA"
    match4_obj_name = "NA"
    match5_obj_name = "NA"

    match1_obj_type = "NA"
    match2_obj_type = "NA"
    match3_obj_type = "NA"
    match4_obj_type = "NA"
    match5_obj_type = "NA"

    try:
        custom_simbad = Simbad()
        custom_simbad.add_votable_fields("otype")
        simbad_result_table = custom_simbad.query_region(
            coord.SkyCoord(ra=my_ra, dec=my_dec, unit=(u.deg, u.deg), frame="fk5"),
            radius=0.83 * u.deg,
        )
        print(simbad_result_table)
        length = len(simbad_result_table)
        loop_lenth = 5
        if length > 5:
            loop_length = 5
        else:
            loop_length = length

        print("length=", length, "loop_length=", loop_length)
        for i in range(loop_length):
            c1 = coord.SkyCoord(ra=my_ra, dec=my_dec, unit=(u.deg, u.deg), frame="fk5")
            RA = simbad_result_table[i]["RA"]
            DEC = simbad_result_table[i]["DEC"]
            OBJ_NAME = simbad_result_table[i]["MAIN_ID"]
            c2 = SkyCoord(RA, DEC, unit=(u.hourangle, u.deg), frame="icrs")
            SEP = c1.separation(c2)
            SEP_ARCSEC = round(SEP.arcsecond, 3)
            OBJ_TYPE = simbad_result_table[i]["OTYPE"]
            print(
                "RA=",
                RA,
                "DEC=",
                DEC,
                "SEP=",
                SEP_ARCSEC,
                " obj_name= ",
                OBJ_NAME,
                "obj_type=",
                OBJ_TYPE,
            )
            if SEP_ARCSEC < 30.01:
                if i == 0:
                    match1_sep = str(SEP_ARCSEC)
                    match1_obj_name = OBJ_NAME
                    match1_obj_type = OBJ_TYPE
                if i == 1:
                    match2_sep = str(SEP_ARCSEC)
                    match2_obj_name = OBJ_NAME
                    match2_obj_type = OBJ_TYPE
                if i == 2:
                    match3_sep = str(SEP_ARCSEC)
                    match3_obj_name = OBJ_NAME
                    match3_obj_type = OBJ_TYPE
                if i == 3:
                    match4_sep = str(SEP_ARCSEC)
                    match4_obj_name = OBJ_NAME
                    match4_obj_type = OBJ_TYPE
                if i == 4:
                    match5_sep = str(SEP_ARCSEC)
                    match5_obj_name = OBJ_NAME
                    match5_obj_type = OBJ_TYPE

    except:
        print("Didn't find match")

    return (
        match1_sep,
        match1_obj_name,
        match1_obj_type,
        match2_sep,
        match2_obj_name,
        match2_obj_type,
        match3_sep,
        match3_obj_name,
        match3_obj_type,
        match4_sep,
        match4_obj_name,
        match4_obj_type,
        match5_sep,
        match5_obj_name,
        match5_obj_type,
    )


# columnn numbers in triple_match_obs.csv
col_no_ra = 14
col_no_dec = 15


with open("triple_match_non_bright.csv", "r") as csvinput:
    with open("triple_match_obs_sim.csv", "w") as csvoutput:
        writer = csv.writer(csvoutput, lineterminator="\n")
        reader = csv.reader(csvinput)
        all = []
        row = next(reader)
        row.append("SIMBAD_match_sep_1")
        row.append("SIMBAD_match_name_1")
        row.append("SIMBAD_match_obj_type_1")
        row.append("SIMBAD_match_sep_2")
        row.append("SIMBAD_match_name_2")
        row.append("SIMBAD_match_obj_type_2")
        row.append("SIMBAD_match_sep_3")
        row.append("SIMBAD_match_name_3")
        row.append("SIMBAD_match_obj_type_3")
        all.append(row)

        for row in reader:
            ra = row[col_no_ra]
            dec = row[col_no_dec]
            (
                match1_sep,
                match1_obj_name,
                match1_obj_type,
                match2_sep,
                match2_obj_name,
                match2_obj_type,
                match3_sep,
                match3_obj_name,
                match3_obj_type,
                match4_sep,
                match4_obj_name,
                match4_obj_type,
                match5_sep,
                match5_obj_name,
                match5_obj_type,
            ) = find_simbad(ra, dec)
            print(match1_sep)
            row.append(match1_sep)
            row.append(match1_obj_name)
            row.append(match1_obj_type)
            row.append(match2_sep)
            row.append(match2_obj_name)
            row.append(match2_obj_type)
            row.append(match3_sep)
            row.append(match3_obj_name)
            row.append(match3_obj_type)
            all.append(row)
        writer.writerows(all)
