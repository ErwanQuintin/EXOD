"""
EXODUS - EPIC XMM-Newton Outburst Detector Ultimate System
Maitrayee Gupta (2022) - maitrayee.gupta@irap.omp.eu

DETECTOR main programme.
"""
import argparse
import os
import sys
import time
from functools import partial
from math import *
from multiprocessing import Pool

import numpy as np

import file_names as FileNames
from file_utils import *
from fits_extractor import *
from logger import logger
from renderer import *
from variability_utils import *

parser = argparse.ArgumentParser()
# Path to files
parser.add_argument("-path", help="Path to the folder containing the observation files", type=str)
parser.add_argument("-out", help="Path to the folder where the output files will be stored", default=None, type=str)

# Variability parameters
parser.add_argument("-bs", "--box-size", dest="bs", help="Size of the detection box in pixel^2.\nDefault: 3", default=3, nargs="?", type=int)
parser.add_argument("-dl", "--detection-level", dest="dl", help="The number of times the median variability is required to trigger a detection.\nDefault: 10", default=10, nargs="?", type=float)
parser.add_argument("-tw", "--time-window", dest="tw", help="The duration of the time windows.\n Default: 100", default=100.0, nargs="?", type=float)
parser.add_argument("-gtr", "--good-time-ratio", dest="gtr", help="Ratio of acceptability for a time window. Shall be between 0.0 and 1.0.\nDefault: 1.0", default=1.0, nargs="?", type=float)
parser.add_argument("-mta", "--max-threads-allowed", dest="mta", help="Maximal number of CPUs the program is allowed to use.\nDefault: 8", nargs="?", default=8, type=int)

# Arguments set by default
parser.add_argument( "-creator", dest="creator", help="User creating the variability files", nargs="?", default=os.environ["USER"], type=str)
parser.add_argument( "-obs", "--observation", dest="obs", help="Observation ID", default=None, nargs="?", type=str)
parser.add_argument( "-inst", "--instrument", dest="inst", help="Type of detector", default="PN", nargs="?", type=str)

# Boolean flags
parser.add_argument("--render", help="Plot variability output, produce pdf", action="store_true")
parser.add_argument("--ds9", help="Plot variability output in emerging ds9 window", action="store_true")
parser.add_argument("--novar", help="Skip variability computation if already done", action="store_true")

args = parser.parse_args()

# Modifying arguments
if args.path[-1] != "/":
    args.path = args.path + "/"
if args.out != None and args.out[-1] != "/":
    args.out = args.out + "/"
if args.out == None:
    args.out = args.path + "{}_{}_{}_{}_{}/".format(
        int(args.dl), int(args.tw), args.bs, args.gtr, args.inst
    )

if args.inst == "PN":
    args.evts = args.path + FileNames.CLEAN_FILE_PN
    args.gti = args.path + FileNames.GTI_FILE_PN
    args.img = args.path + FileNames.IMG_FILE_PN
    ccdnb = 12
if args.inst == "M1":
    args.evts = args.path + FileNames.CLEAN_FILE_M1
    args.gti = args.path + FileNames.GTI_FILE_M1
    args.img = args.path + FileNames.IMG_FILE_M1
    ccdnb = 7
if args.inst == "M2":
    args.evts = args.path + FileNames.CLEAN_FILE_M2
    args.gti = args.path + FileNames.GTI_FILE_M2
    args.img = args.path + FileNames.IMG_FILE_M2
    ccdnb = 7


def main_fct():
    """
    Main function of the detector
    """
    ###
    # Preliminaries
    ###
    print(vars(args))

    print(
        """
        INSTRUMENT      = {0}
        DETECTION LEVEL = {1}
        TIME WINDOW     = {2}
        BOX SIZE        = {3}
        GOOD TIME RATIO = {4}
        """.format(
            args.inst, args.dl, args.tw, args.bs, args.gtr
        )
    )

    print(" Writing output to folder '{0}'".format(args.out))

    # Counter for the overall execution time
    original_time = time.time()

    # Opening the output files
    log_f, var_f, reg_f, best_match_f = open_files(args.out)
    original = sys.stdout
    sys.stdout = Tee(sys.stdout, log_f)

    ###
    # Skipping variability computation if already done
    ###
    vf = False
    if args.novar:
        print(" Checking if variability has been computed.")
        var_f = args.out + FileNames.VARIABILITY
        vf = os.path.isfile(var_f)
        if vf:
            print(" Using existing variability file {0}".format(var_f))
        else:
            print(" No variability file. Applying detector.")
    ###
    # Starting variability computation
    ###

    if not args.novar and not vf:
        # Recovering the EVENTS list
        print(
            " Recovering the events list\t {:7.2f} s".format(
                time.time() - original_time
            )
        )
        try:
            data, header = extraction_photons(args.evts)

            if args.obs == None:
                args.obs = header["OBS_ID"]

        except Exception as e:
            print(" !!!!\nImpossible to extract photons. ABORTING.")
            close_files(log_f, var_f, reg_f, best_match_f)
            exit(-2)

        # Parameters ready
        params = {
            "CREATOR": args.creator,
            "DATE": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
            "OBS_ID": args.obs,
            "INST": args.inst,
            "TW": args.tw,
            "GTR": args.gtr,
            "DL": args.dl,
            "BS": args.bs,
        }

        # Recovering GTI list
        try:
            print(" Extracting data\t\t {:7.2f} s".format(time.time() - original_time))
            gti_list = extraction_deleted_periods(args.gti)

        except Exception as e:
            print(" !!!!\nImpossible to extract gti. ABORTING.")
            close_files(log_f, var_f, reg_f, best_match_f)
            exit(-2)

        # Computation of initial and final time
        t0_observation = min([evt["TIME"] for ccd in data for evt in ccd])
        tf_observation = max([evt["TIME"] for ccd in data for evt in ccd])

        print(
            " Computing variability\t\t {:7.2f} s".format(time.time() - original_time)
        )
        v_matrix = []

        var_calc_partial = partial(
            variability_computation,
            gti_list,
            args.tw,
            args.gtr,
            t0_observation,
            tf_observation,
            args.bs,
            args.inst,
        )

        with Pool(args.mta) as p:
            v_matrix = p.map(var_calc_partial, data)

        # Checking data mode acquisition
        submode = header["SUBMODE"]

        # Applying CCD and Mode configuration
        if args.inst == "PN":
            data_v = PN_config(v_matrix)
            data_v = np.array(data_v)
            if submode == "PrimeLargeWindow":
                data_vm = data_v[:, 100:300]
            elif submode == "PrimeSmallWindow":
                data_vm = data_v[128:192, 200:264]
            else:
                data_vm = data_v

        elif args.inst == "M1":
            data_v = M1_config(v_matrix)
            data_vm = np.array(data_v)

        elif args.inst == "M2":
            data_v = M2_config(v_matrix)
            data_vm = np.array(data_v)

        # Applying geometrical transformation
        if args.inst == "PN":
            img_v = data_transformation_PN(data_vm, header)
        elif args.inst == "M1" or "M2":
            img_v = data_transformation_MOS(data_vm, header)
            # print ("m1 here 2")

        ###
        # Detecting variable areas and sources
        ###

        print(
            " Detecting variable sources\t {:7.2f} s".format(
                time.time() - original_time
            )
        )
        v_matrix = np.array(v_matrix)
        median = np.median(v_matrix)

        # Avoiding a too small median value for detection
        print("\n\tMedian\t\t{0}".format(median))
        if median < 0.75:
            median = 0.75
            print(" \tMedian switched to 0.75. \n")

        variable_areas = []

        # Currying the function for the pool of threads
        variable_areas_detection_partial = partial(
            variable_areas_detection, median, args.bs, args.dl, args.inst
        )
        print("\tBox counts\t{0}".format(args.dl * ((args.bs**2))))
        # Performing parallel detection on each CCD
        with Pool(args.mta) as p:
            variable_areas = p.map(variable_areas_detection_partial, v_matrix)
        sources = variable_sources_position(
            variable_areas,
            args.obs,
            args.inst,
            args.path,
            reg_f,
            log_f,
            args.img,
            best_match_f,
            args.tw,
            args.dl,
            args.bs,
            args.gtr,
        )

        print("\tNb of sources\t{0}\n".format(len(sources)))

        # Writing data to fits file
        fits_writer(img_v, sources, args.img, params, var_f)

    # Plotting variability
    # Renderer
    if args.render:
        print(
            " Rendering variability image\t {:7.2f} s".format(
                time.time() - original_time
            )
        )

        render_variability(
            var_f, args.out + FileNames.OUTPUT_IMAGE, sources=False, maximum_value=None
        )
        render_variability(
            var_f,
            args.out + FileNames.OUTPUT_IMAGE_SRCS,
            sources=True,
            maximum_value=None,
        )

    # ds9
    if args.ds9:
        ds9_renderer(var_f, reg_f)

    print(
        " # Total execution time OBS {0} : {1:.2f} s\n".format(
            args.obs, (time.time() - original_time)
        )
    )
    log_f.close()


if __name__ == "__main__":
    original_time = time.time()
    original_date = time.strftime("%d/%m/%Y %H:%M:%S", time.gmtime())

    main_fct()
