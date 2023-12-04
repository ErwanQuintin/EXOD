"""
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System 
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu       

DETECTOR utilities.

Implementation of variability-related procedures specified into the documentation.
"""
import os
from itertools import combinations
from math import *

import numpy as np
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.simbad import Simbad
from scipy.sparse import csr_matrix
from scipy.stats import linregress

from exodus_utils import check_multiple_sources
from source import Source
from logger import logger


def variability_computation(gti, time_interval, acceptable_ratio, start_time, end_time, inst, box_size, data):
    """
    Calculate variability using the average technique.

    Parameters:
        gti (list): The list of time windows cut-off for the observation.
        time_interval (float): The duration of a time window.
        acceptable_ratio (float): The acceptability ratio for a time window - good time ratio.
        start_time (float): The starting instant of the observation.
        end_time (float): The ending instant of the observation.
        inst (str): Type of detector.
        box_size (float): Box size for the calculation.
        data (list): Events sorted by their TIME attribute.

    Returns:
        list: The matrix V_round.
    """
    logger.debug(f'gti: {gti}')
    logger.debug(f'time_interval: {time_interval}')
    logger.debug(f'acceptable_ratio: {acceptable_ratio}')
    logger.debug(f'start_time: {start_time}')
    logger.debug(f'end_time:{end_time}')
    logger.debug(f'inst:{inst}')
    logger.debug(f'box_size:{box_size}')
    logger.debug('data:')
    logger.debug(data)

    # Defining the variables and matrices
    if inst == "PN":
        xMax = 64
        yMax = 200
    elif inst == "M1" or "M2":
        xMax = 600
        yMax = 600
    
    n_bins = int(np.ceil((end_time - start_time) / time_interval))
    stop_time = start_time + n_bins * time_interval
    time_ratio = (stop_time - end_time) / time_interval 

    logger.debug(f'n_bins={n_bins} stop_time={stop_time} time_ratio={time_ratio:.2f}')
    if time_ratio > acceptable_ratio:
        logger.debug(f'time ratio ({time_ratio:.2f} > acceptable_ratio {acceptable_ratio:.2f}')
        n_bins = n_bins - 1
        stop_time = start_time + n_bins * time_interval
        logger.debug('Reducing nbins by 1. new stop_time={stop_time}')

    time_windows = np.arange(start_time, stop_time, time_interval)
    logger.debug(f'number of time_windows={len(time_windows)}')

    # We treat the GTIs, depending on whether we accept partial time bins or not
    if acceptable_ratio < 1:
        logger.debug('acceptable_ratio<1, only accepting partial time bins')
        ##If we accept partial time bins, the effective fraction of GTI in a time bin is computed by oversampling in each
        ##time bin, and getting the fraction of these subsamples that are in GTIs
        oversampling = 1000
        oversampled_timewindows = np.arange(start_time, stop_time, time_interval / oversampling)

        indices_oversampledtimebins_gtistart = np.searchsorted(gti["START"], oversampled_timewindows[:-1])
        indices_oversampledtimebins_gtistop  = np.searchsorted(gti["STOP"], oversampled_timewindows[1:])

        indices_fully_inGTI = np.where(indices_oversampledtimebins_gtistart - indices_oversampledtimebins_gtistop==1,1,0)
        splitted_oversample = [indices_fully_inGTI[i : i + oversampling] for i in range(0, len(indices_fully_inGTI), oversampling)]
        projection_ratio    = np.array([np.sum(elt) / len(elt) for elt in splitted_oversample])

    else:
        logger.debug('acceptable_ratio=1, only allowing full time bins')
        ##If we only accept full GTI time bins, we check if the beginning and end of each time bin are in the same GTI,
        ##by looking at the difference of insertion index. If it is 1, then the GTI is the same

        indices_timebins_gtistart = np.searchsorted(gti["START"], time_windows[:-1])
        indices_timebins_gtistop = np.searchsorted(gti["STOP"], time_windows[1:])
        projection_ratio = np.where(indices_timebins_gtistart - indices_timebins_gtistop == 1, 1, np.nan)

        logger.debug('indices_timebins_gtistart:')
        logger.debug(indices_timebins_gtistart)
        logger.debug('indices_timebins_gtistop:')
        logger.debug(indices_timebins_gtistop)
        logger.debug('projection_ratio:')
        logger.debug(projection_ratio)

    # Counting events
    ## Group by pixel and time, first on X then on Y and t, selecting boxes at the same time
    nbr_events = len(data)
    box_halfwidth = (box_size - 1) // 2  # Only works well for odd box size. If even, box_size-1 is used
    grouped_events = [data[np.where(np.abs(data["RAWX"] - x) <= box_halfwidth)] for x in range(xMax)]
    grouped_events = [[dataX[np.where(np.abs(dataX["RAWY"] - y) <= box_halfwidth)] for y in range(yMax)] for dataX in grouped_events]
    grouped_events = [[[evt["TIME"] for evt in grouped_events[x][y]] for y in range(yMax)] for x in range(xMax)]

    logger.debug(f'nbr_events: {nbr_events} box_halfwidth={box_halfwidth} len(grouped_events)={len(grouped_events)}')

    ## If the matrix is sparse (i.e. less than half dense) we use sparse matrices, otherwise remain arrays
    cdt = np.where(projection_ratio >= acceptable_ratio)[0]
    logger.debug('cdt:')
    logger.debug(cdt)

    if nbr_events < 0.5 * (xMax * yMax * len(time_windows)):
        # This creates a 2D matix of 1D sparse (time) vectors, built by np.histogram on the photon times
        # This takes a long time
        logger.debug('Getting Counted events')
        counted_events = np.array([[csr_matrix(np.histogram(grouped_events[x][y], time_windows)[0][cdt] / projection_ratio[cdt]) for y in range(yMax)] for x in range(xMax)])
        logger.debug('Calculating Maxes')
        image_max    = np.array([[c.max() for c in column] for column in counted_events])
        logger.debug('Calculating minimums')
        image_min    = np.array([[c.min() for c in column] for column in counted_events])
        logger.debug('Calculating medians')
        image_median = np.array([[c.mean() for c in column] for column in counted_events])


        #logger.debug('counted_events:')
        #logger.debug(counted_events)

        logger.debug('image_max:')
        logger.debug(image_max)

        logger.debug('image_min:')
        logger.debug(image_min)

        logger.debug('image_median:')
        logger.debug(image_median)

    else:
        ##If the matrix is not sparse, we create this 3D (XxYxTime) array
        counted_events = np.array([[np.histogram(grouped_events[x][y], time_windows)[0][cdt] / projection_ratio[cdt] for y in range(yMax)] for x in range(xMax)], dtype=int)
        image_max = np.max(counted_events, axis=2)
        image_min = np.min(counted_events, axis=2)
        image_median = np.median(counted_events, axis=2)

    # Computing variability
    V_mat = np.where(image_median > 0, np.max((image_max - image_median, image_median - image_min)) / image_median, image_max)
    return V_mat


########################################################################
#                                                                      #
#            Detecting variable areas                                  #
#                                                                      #
########################################################################


def box_computations(variability_matrix, x, y, box_size, inst):
    """
    Sum the variability values within a specified box.

    Parameters:
    - variability_matrix (list of lists): The V round matrix.
    - x (int): The x-coordinate of the top-left corner of the box.
    - y (int): The y-coordinate of the top-left corner of the box.
    - box_size (int): The length of a side of the box.
    - inst (str): The instrument type ("PN", "M1", or "M2").

    Returns:
    int: The sum of the variability for each pixel within the box.

    Raises:
    AssertionError: If the box is out of the limits of the CCD.
    """
    if inst == "PN":
        assert x <= 63 - box_size
        assert y <= 199 - box_size
    elif inst == "M1" or "M2":
        assert x <= 599 - box_size
        assert y <= 599 - box_size
    # Exception raised if box out of the limits of the CCD

    cpt = 0
    for i in range(x, x + box_size):
        for j in range(y, y + box_size):
            cpt += variability_matrix[x][y]

    return cpt


########################################################################
#                                                                      #
#  Function summing the variability values of a box                    #
#                                                                      #
########################################################################


def __add_to_detected_areas(x, y, box_size, detected_areas, variability_matrix):
    """
    Function summing the variability values into a box.
    @param x:   The x coordinate of the top-left corner of the box
    @param y:   The y coordinate of the top-left corner of the box
    @param box_size:   The length of a side of a box
    @param detected_areas:  The A round set, containing already detected areas
    @return: The sum of the variability for each pixel of the box ??? Should be the A round set, containing the updated detected areas
    """
    box_set = {
        (a, b, variability_matrix[a][b])
        for a in range(x, x + box_size)
        for b in range(y, y + box_size)
    }
    inserted = False
    i = 0
    while (not inserted) and i < len(detected_areas):
        # If a part of the box has already been detected, the two sets of coordinates are merged...
        if len(box_set & detected_areas[i]) > 1:
            detected_areas[i] |= box_set
            # Equivalent to detected_areas[i] = detected_areas[i] | box_set => add area of the box to the detected areas
            # ...and the loop is over
            inserted = True
        i += 1

    # If there has not been any merger :
    if not inserted:
        detected_areas.append(box_set)
    return detected_areas


########################################################################
#                                                                      #
#   Remove objects that are part of a readout streak                   #
#                                                                      #
########################################################################


def remove_readout_streak(input_table):
    """
    Function to emove objects that are part of a readout streak
    @param input_table:  Input trable with RA,DEC positions
    @return: An outout table with t objects that are part of a readout streak removed.
    """
    cut_off_err = 0.1
    cut_off_dist = 0.008

    ra = input_table["RA"][:]
    dec = input_table["DEC"][:]
    if len(ra) < 6:
        return input_table
    xy = []
    for i in range(0, len(ra)):
        xy.append([ra[i], dec[i]])

    comb = combinations(xy, 6)
    fit_slope = -100
    fit_intercept = -100
    x_new = [min(ra), max(ra)]
    y_new = []
    for i in list(comb):
        x_arr = []
        y_arr = []
        for j in range(0, len(i)):
            x_arr.append(i[j][0])
            y_arr.append(i[j][1])
        slope, intercept, r_value, p_value, std_err = linregress(x_arr, y_arr)
        if std_err < cut_off_err:
            fit_slope = slope
            fit_intercept = intercept
            y_new.append((slope * x_new[0]) + intercept)
            y_new.append((slope * x_new[1]) + intercept)
            break
    logger.debug(f'fit_slope={fit_slope}')
    logger.debug(f'fit_intercept={fit_intercept}')
    if (fit_slope != -100) and (fit_intercept != -100):
        logger.debug("found hit")

    remove_ids = []
    remove_ra = []
    remove_dec = []
    for it in range(0, len(ra)):
        p1 = np.array([x_new[0], y_new[0]])
        p2 = np.array([x_new[1], y_new[1]])
        p3 = np.array([ra[it], dec[it]])
        d = abs(np.cross(p2 - p1, p3 - p1) / np.linalg.norm(p2 - p1))
        if d < cut_off_dist:
            remove_ids.append(it)
            remove_ra.append(ra[it])
            remove_dec.append(dec[it])

    result_table = Table(
        names=(
            "ID",
            "INST",
            "CCDNR",
            "RAWX",
            "RAWY",
            "RAWR",
            "X",
            "Y",
            "SKYR",
            "RA",
            "DEC",
            "R",
            "VCOUNT",
        ),
        dtype=(
            "i2",
            "U25",
            "i2",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
        ),
    )

    for it in range(0, len(input_table)):
        found = 0
        for it2 in range(0, len(remove_ra)):
            if (input_table["RA"][it] == remove_ra[it2]) and (
                input_table["DEC"][it] == remove_dec[it2]
            ):
                found = found + 1
        if found == 0:
            result_table.add_row(input_table[it])

    for it in range(0, len(result_table)):
        result_table["ID"][it] = it

    return result_table



def variable_areas_detection(lower_limit, box_size, detection_level, inst, variability_matrix):
    """
    Detects variable areas in a variability matrix.

    Args:
        lower_limit (float): The smallest variability value needed to consider a pixel variable.
        box_size (int, optional): The size of the box (default is 3).
        detection_level (float): A factor for the limit of detection.
        inst: Unused parameter, reserved for future use.
        variability_matrix (list of lists): The matrix returned by variability_calculation.

    Returns:
        list of sets: A list of sets of coordinates for each area detected as variable.
    """
    logger.debug(f'lower_limit={lower_limit} box_size={box_size} detection_level={detection_level} inst={inst}')
    
    output = []
    box_count = 0

    for i in range(len(variability_matrix) - box_size):
        j = 0
        m = 0  # boxes above detection level counter
        while j < len(variability_matrix[i]) - box_size:
            box_count = box_computations(variability_matrix, i, j, box_size, inst)

            # If there's nothing into the box, it is completely skipped
            if box_count == 0:
                j += box_size
            else:
                box_count_threshold = detection_level * ((box_size**2) * lower_limit)
                logger.debug(f'i={i} j={j} box_count:{box_count} box_count_threshold={box_count_threshold}')            
                if box_count > box_count_threshold:
                    output = __add_to_detected_areas(i, j, box_size, output, variability_matrix)
                    m += 1
                j += 1
    return output

def variable_sources_position(variable_areas_matrix, obs, inst, path_out, reg_file, img_file, best_match_file, bs, dl, tw, gtr):
    """
    Compute the position of the detected variable sources.

    Parameters:
        variable_areas_matrix (numpy.ndarray): A matrix representing the detected variable areas.
        obs (str): The EPIC-pn OBSID, which will be written in the output file.
        inst (str): Information about the instrument used.
        path_out (str): The path to the output directory where result files will be saved.
        reg_file (str): The file path for the region file where the sources will be written.
        img_file (str): The file path for the image file.
        best_match_file (str): The file path for the best match file.
        bs (float): The size of the box used in the computation.
        dl (float): The Detection level value used in the computation.
        tw (float): The time window value used in the computation.
        gtr (float): The good time ratio value used in the computation.

    Returns:
        astropy.table.Table: A table object containing the parameters of the detected sources.
    """
    sources = []
    cpt_source = 0
    correction_factor_x = 3
    correction_factor_y = 3
    region_number = 1

    #TODO looks like the correction factors have been set to 0, ??
    if inst == "PN":
        ccdnb = 12
        correction_factor_x = 0  # 3/2
        correction_factor_y = 0  # 3/2
    elif inst == "M1" or "M2":
        ccdnb = 7
        correction_factor_x = 0  # 3
        correction_factor_y = 0  # 3

    logger.debug('Computing Source position...')
    for ccd in range(ccdnb):
        for source in variable_areas_matrix[ccd]:
            center_x = 0
            center_y = 0
            it = 0
            max_ = 0
            for p in source:
                if p[2] > max_:
                    center_x = p[0]
                    center_y = p[1]
                    max_ = p[2]
                it = it + 1
            logger.debug(f'it={it} center_x={center_x} center_y={center_y} max_={max_}')
            # Sweet jesus...
            r = round(sqrt((max([abs((p[0] + correction_factor_x) - center_x) for p in source]))**2
                    + (max([abs((p[1] + correction_factor_y) - center_y) for p in source]))**2),2)
            vcount = max_

            # Avoiding bad pixels
            bad_pixels = [["PN", 4, 11],
                          ["PN", 4, 12],
                          ["PN", 4, 13],
                          ["PN", 5, 12],
                          ["PN", 10, 28],
                          ["M1", 1, 318]]
            
            logger.debug('Avoiding bad pixels...')
            if [inst, ccd + 1, int(center_x)] not in bad_pixels:
                cpt_source += 1
                sources.append([cpt_source, inst, ccd + 1, center_x, center_y, r, vcount])
                source = Source(cpt_source, inst, ccd+1, center_x, center_y, r)
    

    logger.debug('Creating source table...')
    source_table = Table(
        names=(
            "ID",
            "INST",
            "CCDNR",
            "RAWX",
            "RAWY",
            "RAWR",
            "X",
            "Y",
            "SKYR",
            "RA",
            "DEC",
            "R",
            "VCOUNT",
        ),
        dtype=(
            "i2",
            "U25",
            "i2",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
            "f8",
        ),
    )

    logger.debug("Filling the source table...")
    for s in sources:
        src = Source(id_src=s[0], inst=s[1], ccd=s[2]+1, rawx=s[3], rawy=s[4], rawr=s[5])
        src.sky_coord(path_out, img_file)

        src.r = round(src.r, 3)
        # Adding source to table
        source_table.add_row(
            [
                src.id_src,
                src.inst,
                src.ccd,
                src.rawx,
                src.rawy,
                src.rawr,
                src.x,
                src.y,
                src.skyr,
                src.ra,
                src.dec,
                src.r,
                src.vcount,
            ]
        )

    # Removing multiple sources
    source_table = check_multiple_sources(source_table)
    source_table = remove_readout_streak(source_table)



#    logger.debug(f'Writing header to {best_match_file}')
#    with open(best_match_file, "w+") as f:
#        f.write("OBS_ID,REGION_NUMBER,SIMBAD_MATCH_NO,OBJ_NAME,RA_EXOD,DEC_EXOD,RA_SIM,DEC_SIM,SEP,OBJ_TYPE,INST,DL,GTR,BS,TW\n")
#    
#    logger.debug(f'Writing header to Master_Catalogue.csv')
#    with open('Master_Catalogue.csv', 'w+') as f:
#        f.write("OBS_ID,REGION_NUMBER,SIMBAD_MATCH_NO,OBJ_NAME,RA_EXOD,DEC_EXOD,RA_SIM,DEC_SIM,SEP,OBJ_TYPE,INST,DL,GTR,BS,TW\n")
#
#
#
#    # Head text
#    text = """# Region file format: DS9 version 4.0 global
#    # XMM-Newton OBSID {0}
#    # Instrument {1}
#    # EXOD variable sources
#    
#    global color=green font="times 8 normal roman"
#    j2000
#    
#    """.format(obs, inst)
#
#    # ds9 text
#    for s in source_table:
#        logger.debug(s)
#
#
#
#    for s in source_table:
#        text = text + 'circle {0}, {1}, {2}" # text="{3}"\n'.format(s["RA"], s["DEC"], s["R"], s["ID"])
#        try:
#            custom_simbad = Simbad()
#            custom_simbad.add_votable_fields("otype")
#            simbad_result_table = custom_simbad.query_region(
#                coord.SkyCoord(
#                    ra=s["RA"], dec=s["DEC"], unit=(u.deg, u.deg), frame="fk5"
#                ),
#                radius=0.0083 * u.deg,
#            )
#
#            length = len(simbad_result_table)
#            loop_lenth = 5
#            if length > 5:
#                loop_length = 5
#            else:
#                loop_length = length
#
#            for i in range(loop_length):
#                c1 = coord.SkyCoord(
#                    ra=s["RA"], dec=s["DEC"], unit=(u.deg, u.deg), frame="fk5"
#                )
#                RA_EXOD = s["RA"]
#                DEC_EXOD = s["DEC"]
#                RA = simbad_result_table[i]["RA"]
#                DEC = simbad_result_table[i]["DEC"]
#                OBJ_NAME = simbad_result_table[i]["MAIN_ID"]
#                c2 = SkyCoord(RA, DEC, unit=(u.hourangle, u.deg), frame="icrs")
#                SEP = c1.separation(c2)
#                SEP_ARCSEC = round(SEP.arcsecond, 3)
#                OBJ_TYPE = simbad_result_table[i]["OTYPE"]
#
#                line = """{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n""".format(
#                    obs,
#                    region_number,
#                    i + 1,
#                    OBJ_NAME,
#                    RA_EXOD,
#                    DEC_EXOD,
#                    RA,
#                    DEC,
#                    SEP_ARCSEC,
#                    OBJ_TYPE,
#                    inst,
#                    dl,
#                    gtr,
#                    bs,
#                    tw,
#                )
#
#                with best_match_file
#                match_file.write(line)
#                master_file.write(line)
#
#        except:
#            pass
#            line = """{0},{1},N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A\n""".format(
#                obs, region_number
#            )
#            match_file.write(line)
#            master_file.write(line)
#            logger.debug(f"Didn't find match for RA=", {s["RA"]}, "DEC=", {s["DEC"]})
#
#        region_number = region_number + 1
#
#    match_file.close()
#    master_file.close()
#
#    # Writing region file
#    reg_f = open(reg_file, "w")
#    reg_f.write(text)
#    reg_f.close()
    return source_table
